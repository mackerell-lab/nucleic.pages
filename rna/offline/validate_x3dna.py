#!/usr/bin/env python3
"""Independent RNA geometry oracle using a separately installed 3DNA executable.

The fixture coordinate adapter is deliberately limited to canonical, single-altloc
PDB fixtures. Production mmCIF normalization is not validated by this adapter.
Raw downloaded coordinates and 3DNA outputs remain outside the Pages repository.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import urllib.request

RNA_PARAMS = ['chi','alpha','beta','gamma','delta','epsilon','zeta','e_z',
              'eta','theta','eta1','theta1','eta2','theta2',
              'v0','v1','v2','v3','v4','tm','p','sszp','dp']
LABEL = re.compile(r'^\s*\d+\s+([^:]+):([.\w-]+)_:\[([^]]+)\][A-Za-z]')

def nt_key(chain, comp, seq):
    # Numerical fixture ordering prevents lexicographic comp/author-ID ordering
    # from reversing same-chain hairpin pairs in the production adapter.
    return f'{chain.strip()}.{int(seq):06d}.{comp.lstrip(".")}'

def number(text):
    try:
        value = float(text)
        return value if math.isfinite(value) else None
    except ValueError:
        return None

def parse_torsions(path):
    """Parse fixed columns: a blank chi classification must not shift alpha."""
    rows = defaultdict(dict)
    section = None
    for line in path.read_text().splitlines():
        if line.startswith('Main chain and chi'): section = 'backbone'
        elif line.startswith('Pseudo (virtual)'): section = 'pseudo'
        elif line.startswith('Sugar conformational'): section = 'sugar'
        match = LABEL.match(line)
        if not match or not section: continue
        chain, resid, comp = match.groups()
        key = nt_key(chain, comp, resid.lstrip('.'))
        tail = line[match.end():]
        if section == 'backbone':
            rows[key]['chi'] = number(tail[:8])
            # check_chi appends a six-character classification after chi.
            offset = 14
            for i, param in enumerate(['alpha','beta','gamma','delta','epsilon','zeta','e_z']):
                rows[key][param] = number(tail[offset+8*i:offset+8*(i+1)])
        elif section == 'pseudo':
            for i,param in enumerate(['eta','theta','eta1','theta1','eta2','theta2']):
                rows[key][param] = number(tail[8*i:8*(i+1)])
        else:
            for i,param in enumerate(['v0','v1','v2','v3','v4','tm','p']):
                rows[key][param] = number(tail[8*i:8*(i+1)])
            rows[key]['sszp'] = number(tail[-16:-8])
            rows[key]['dp'] = number(tail[-8:])
    if not rows: raise ValueError(f'No torsion rows parsed: {path}')
    return dict(rows)

def fixture_entry(pdb_id, path):
    residues = {}
    model = None
    for line in path.read_text().splitlines():
        if line.startswith('MODEL') and model is None: model = line[10:14].strip()
        if line.startswith('ENDMDL'): break
        if line[:6] not in ('ATOM  ', 'HETATM'): continue
        comp = line[17:20].strip()
        if comp not in {'A','C','G','U'}: continue
        if line[16].strip(): raise ValueError('This oracle adapter does not handle altlocs')
        chain, seq, insertion = line[21].strip() or '_', line[22:26].strip(),line[26].strip()
        if insertion: raise ValueError('Insertion-code fixtures require mmCIF identity mapping')
        key = nt_key(chain, comp, seq)
        if key not in residues:
            residues[key] = {'id':key,'pdb_id':pdb_id,'comp_id':comp,'label_asym_id':chain,
                             'chain_id':chain,'label_seq_id':seq,'auth_seq_id':seq,
                             'model_id':model or '1','atoms':{}}
        atom = line[12:16].strip().replace('*',"'")
        if atom in residues[key]['atoms']: raise ValueError('Duplicate fixture atom')
        residues[key]['atoms'][atom] = [float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]]
    values = list(residues.values()); links=[]
    for left,right in zip(values,values[1:]):
        if left['chain_id'] != right['chain_id']: continue
        if int(right['auth_seq_id']) != int(left['auth_seq_id'])+1: continue
        a,b=left['atoms'].get("O3'"),right['atoms'].get('P')
        if a and b and 0.8 <= math.dist(a,b) <= 2.4:
            links.append({'from_id':left['id'],'to_id':right['id'],'status':'connected'})
    return {'pdb_id':pdb_id,'residues':values,'links':links}

PAIR_PARAMETERS = ['shear','stretch','stagger','buckle','propeller','opening']
QUALITY_PARAMETERS = ['lambda_1','lambda_2','c1c1','rn9_yn1','rc8_yc6']
STEP_PARAMETERS = ['shift','slide','rise','tilt','roll','twist']
HELICAL_PARAMETERS = ['x_disp','y_disp','h_rise','inclination','tip','h_twist']
POSITION_PARAMETERS = ['xp','yp','zp','xph','yph','zph']
RADIUS_PARAMETERS = ['strand_i_p_radius','strand_i_o4_radius','strand_i_c1_radius','strand_ii_p_radius','strand_ii_o4_radius','strand_ii_c1_radius']
SAME_PARAMETERS = ['strand_i_p_p','strand_i_c1_c1','strand_ii_p_p','strand_ii_c1_c1']

def parse_geometry(folder):
    left=re.compile(r'>([^:]+):([.\w-]+)_:\[([^]]+)\][A-Za-z]')
    right=re.compile(r'[A-Za-z]\[([^]]+)\]:([.\w-]+)_:([^<]+)<')
    pairs=[]
    for line in (folder/'pairs.inp').read_text().splitlines():
        a,b=left.search(line),right.search(line)
        if not a or not b: continue
        pairs.append([nt_key(a[1],a[3],a[2].lstrip('.')),nt_key(b[3],b[1],b[2].lstrip('.'))])
    tables={name:{} for name in ['base_pairs','lambda','steps','helical','step_position','same_strand','helix_radius']}
    names={'base_pairs':PAIR_PARAMETERS,'lambda':QUALITY_PARAMETERS,'steps':STEP_PARAMETERS,
           'helical':HELICAL_PARAMETERS,'step_position':POSITION_PARAMETERS,'same_strand':SAME_PARAMETERS,'helix_radius':RADIUS_PARAMETERS}
    section=None
    for line in (folder/'selected.out').read_text().splitlines():
        if line.startswith('Local base-pair parameters'):section='base_pairs'
        elif line.startswith('Local base-pair step parameters'):section='steps'
        elif line.startswith('Local base-pair helical parameters'):section='helical'
        elif line.strip().startswith('bp     lambda(I)'):section='lambda'
        elif line.strip().startswith('step       Xp      Yp      Zp'):section='step_position'
        elif line.startswith('Same strand P--P'):section='same_strand'
        elif line.startswith('Helix radius (radial'):section='helix_radius'
        elif line.startswith(('*','-','Simple ')) or '~' in line:section=None
        if not section or not re.match(r'^\s*\d+\s+',line):continue
        tokens=line.split();index=int(tokens[0])-1
        if index>=len(pairs):raise ValueError('Pair index out of range')
        if section in ['base_pairs','lambda']:key='|'.join(pairs[index])
        else:
            if index+1>=len(pairs):continue
            key='|'.join(pairs[index]+pairs[index+1])
        values=tokens[2:4]+tokens[6:8] if section=='same_strand' else tokens[2:2+len(names[section])]
        if len(values)!=len(names[section]):raise ValueError(f'Unexpected geometry row: {line}')
        tables[section][key]=dict(zip(names[section],map(number,values)))
    return pairs,tables

def analytical_local_reference(entry):
    def subtract(a,b):return [x-y for x,y in zip(a,b)]
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    def norm(a):return math.sqrt(dot(a,a))
    def angle(a,b,c):
        x,y=subtract(a,b),subtract(c,b)
        return math.degrees(math.acos(max(-1,min(1,dot(x,y)/(norm(x)*norm(y))))))
    def dihedral(a,b,c,d):
        x,y,z=subtract(b,a),subtract(c,b),subtract(d,c)
        n1,n2=cross(x,y),cross(y,z)
        return math.degrees(math.atan2(dot(cross(n1,n2),y)/norm(y),dot(n1,n2)))
    result={}
    for residue in entry['residues']:
        purine=residue['comp_id'] in 'AG';n='N9' if purine else 'N1';atoms=residue['atoms']
        definitions={'o4_c1_n':["O4'","C1'",n],'c2_c1_n':["C2'","C1'",n],
                     'c2_o2_length':["C2'","O2'"], 'c1_c2_o2':["C1'","C2'","O2'"],
                     'c3_c2_o2':["C3'","C2'","O2'"],'o4_c1_c2_o2':["O4'","C1'","C2'","O2'"]}
        if purine:definitions.update({'c1_n9_c4':["C1'",n,'C4'],'c1_n9_c8':["C1'",n,'C8']})
        else:definitions.update({'c1_n1_c2':["C1'",n,'C2'],'c1_n1_c6':["C1'",n,'C6']})
        result[residue['id']]={}
        for parameter,definition in definitions.items():
            points=[atoms.get(name) for name in definition]
            value=None
            if all(point is not None for point in points):
                value=math.dist(*points) if len(points)==2 else angle(*points) if len(points)==3 else dihedral(*points)
            result[residue['id']][parameter]=value
    return result

def execute(cmd, cwd, env):
    result=subprocess.run(cmd,cwd=cwd,env=env,text=True,capture_output=True)
    with (cwd/'commands.log').open('a') as file:
        file.write(json.dumps(cmd)+'\n'+result.stdout+result.stderr+'\n')
    if result.returncode: raise RuntimeError(f'Command failed: {cmd}; see {cwd}/commands.log')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--x3dna-root',required=True,type=Path)
    parser.add_argument('--work-dir',required=True,type=Path)
    parser.add_argument('--pdb-id',action='append',default=[])
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args(); work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
    xroot=args.x3dna_root.resolve();env=os.environ.copy();env['X3DNA']=str(xroot);env['PATH']=str(xroot/'bin')+':'+env.get('PATH','')
    entries=[]; refs={}; provenance=[]; geometry_cases=[]; geometry_refs={}
    for pdb_id in (args.pdb_id or ['1rna','1sdr','1t4x','433d','1zih']):
        pdb_id=pdb_id.lower()
        if not re.fullmatch('[0-9a-z]{4}',pdb_id):raise ValueError('Expected a four-character fixture accession')
        folder=work/pdb_id;folder.mkdir(exist_ok=True);raw=folder/'download.pdb'
        url=f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb'
        if not raw.exists():raw.write_bytes(urllib.request.urlopen(url,timeout=60).read())
        selected=[]
        for line in raw.read_text().splitlines():
            if line.startswith('ENDMDL'):break
            if line.startswith(('ATOM  ','HETATM','TER   ')):selected.append(line)
        selected_path=folder/'selected.pdb';selected_path.write_text('\n'.join(selected)+'\nEND\n')
        for command in [['find_pair','selected.pdb','pairs.inp'],['analyze','pairs.inp'],['analyze','-t=residues.tor','selected.pdb']]:
            execute(command,folder,env)
        entries.append(fixture_entry(pdb_id,raw));refs[pdb_id]=parse_torsions(folder/'residues.tor')
        pair_order,geometry_refs[pdb_id]=parse_geometry(folder)
        comps={residue['id']:residue['comp_id'] for residue in entries[-1]['residues']}
        geometry_cases.append({'entry':entries[-1],'graph':{'status':'available','edges':[
            {'id':f'{pdb_id}:pair:{i}','category':'basepair','family':'cWW' if comps[a]+comps[b] in {'AU','UA','GC','CG','GU','UG'} else 'unclassified','residue1_id':a,'residue2_id':b,'near':False,'alternative':False}
            for i,(a,b) in enumerate(pair_order)]}})
        provenance.append({'pdb_id':pdb_id,'source_url':url,'raw_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
                           'selected_sha256':hashlib.sha256(selected_path.read_bytes()).hexdigest()})
    fixtures=work/'entries.json'; fixtures.write_text(json.dumps(entries,indent=2)+'\n')
    (work/'reference.json').write_text(json.dumps(refs,indent=2)+'\n')
    if args.prepare_only:
        print(json.dumps({'prepared':len(entries),'entries':str(fixtures)}));return
    runner=Path(__file__).resolve().parents[1]/'tests/reference/compute_residues.mjs'
    output=subprocess.check_output(['node',str(runner),str(fixtures)],text=True)
    computed=json.loads(output);checks=[]
    for entry in computed:
        for row in entry['rows']:
            reference=refs[entry['pdb_id']].get(row['id'])
            if reference is None: raise ValueError(f'Missing reference residue: {row["id"]}')
            for parameter in RNA_PARAMS:
                actual=row['values'].get(parameter);expected=reference.get(parameter)
                # Explorer defines e-z as the shortest signed difference; 3DNA prints
                # the unwrapped difference of epsilon360 and zeta360.
                if parameter == 'e_z' and expected is not None: expected=(expected+180)%360-180
                tolerance=0.006 if parameter in {'sszp','dp'} else 0.051
                error=None
                if actual is None and expected is None: status='both_missing'
                elif actual is None or expected is None: status='availability_mismatch'
                else:
                    error=abs(actual-expected)
                    if parameter not in {'tm','sszp','dp','e_z'}:error=abs((actual-expected+180)%360-180)
                    status='pass' if error<=tolerance else 'fail'
                checks.append({'pdb_id':entry['pdb_id'],'id':row['id'],'parameter':parameter,
                               'actual':actual,'reference':expected,'error':error,'tolerance':tolerance,'status':status})
    analytical={entry['pdb_id']:analytical_local_reference(entry) for entry in entries}
    for entry in computed:
        for row in entry['rows']:
            for parameter,expected in analytical[entry['pdb_id']][row['id']].items():
                actual=row['values'].get(parameter);error=None if actual is None or expected is None else abs(actual-expected)
                status='both_missing' if actual is None and expected is None else 'pass' if error is not None and error<1e-8 else 'fail'
                checks.append({'pdb_id':entry['pdb_id'],'id':row['id'],'parameter':parameter,'oracle':'independent_python_atom_formula',
                               'actual':actual,'reference':expected,'error':error,'tolerance':1e-8,'status':status})
    geometry_input=work/'geometry_inputs.json';geometry_input.write_text(json.dumps(geometry_cases)+'\n')
    geometry_runner=runner.with_name('compute_geometry.mjs')
    geometry=json.loads(subprocess.check_output(['node',str(geometry_runner),str(geometry_input)],text=True))
    (work/'geometry_reference.json').write_text(json.dumps(geometry_refs,indent=2)+'\n')
    (work/'geometry_actual.json').write_text(json.dumps(geometry,indent=2)+'\n')
    for entry in geometry:
        for level,groups in [('pairs',['base_pairs','lambda']),('steps',['steps','helical','step_position','same_strand','helix_radius'])]:
            for row in entry[level]:
                key='|'.join(row['residue_ids'])
                for group in groups:
                    reference=geometry_refs[entry['pdb_id']][group].get(key)
                    if reference is None:raise ValueError(f'Missing {group} reference for {key}')
                    for parameter,expected in reference.items():
                        actual=row['values'].get(parameter);tolerance=.051 if group=='lambda' else .006
                        error=None if actual is None or expected is None else abs(actual-expected)
                        if error is not None and parameter in ['buckle','propeller','opening','tilt','roll','twist','inclination','tip','h_twist']:
                            error=abs((actual-expected+180)%360-180)
                        status='both_missing' if actual is None and expected is None else 'pass' if error is not None and error<=tolerance else 'fail'
                        checks.append({'pdb_id':entry['pdb_id'],'id':key,'parameter':parameter,'group':group,'oracle':'x3dna_analyze',
                                       'actual':actual,'reference':expected,'error':error,'tolerance':tolerance,'status':status})
    counts=dict(Counter(check['status'] for check in checks))
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'fixture_count':len(entries),
            'residue_count':sum(len(x['residues']) for x in entries),'counts':counts,'parameters':sorted({check['parameter'] for check in checks}),
            'provenance':provenance,'x3dna_root':str(xroot),
            'x3dna_analyze_sha256':hashlib.sha256((xroot/'bin/analyze').read_bytes()).hexdigest(),
            'reference_source_sha256':{str(path.relative_to(xroot)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in [xroot/'src/ana_fncs.c',xroot/'src/cmn_fncs.c',*[xroot/f'config/Atomic_{base}.pdb' for base in 'ACGU']]},
            'scope':'Canonical first-model PDB fixture residue/pair/step geometry and independent RNA atom formulas. Pair identities supplied by 3DNA, not the production interaction detector. Not production mmCIF normalization, all-PDB discovery, FR3D detection, modified RNA, or mixed-helicity stem validation.',
            'checks':checks}
    (work/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
    if counts.get('fail',0)+counts.get('availability_mismatch',0):raise SystemExit(1)

if __name__=='__main__':main()
