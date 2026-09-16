import { readFile } from 'node:fs/promises';
import { computeGeometry } from '../../offline/geometry_adapter.mjs';
const fixtures=JSON.parse(await readFile(process.argv[2],'utf8'));
process.stdout.write(JSON.stringify(fixtures.map(({entry,graph})=>({pdb_id:entry.pdb_id,...computeGeometry(entry,graph)}))));
