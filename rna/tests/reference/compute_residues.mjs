import { readFile } from 'node:fs/promises';
import { computeResidueObservables } from '../../offline/residue_geometry.mjs';
const entries=JSON.parse(await readFile(process.argv[2], 'utf8'));
const output=entries.map(entry=>({pdb_id:entry.pdb_id,rows:computeResidueObservables(entry)}));
process.stdout.write(JSON.stringify(output));
