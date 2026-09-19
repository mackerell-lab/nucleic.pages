import { RnaDataRepository } from './core/repository.js';
import { PureRnaExplorer } from './app/PureRnaExplorer.js';
import { mountStaticControlHelp } from './views/panels.js';

const root = document.querySelector('#rnaExplorer');
mountStaticControlHelp(root);
const repository = new RnaDataRepository({ manifestUrl: new URL('../assets/pure_rna/manifest.json', import.meta.url) });
const explorer = new PureRnaExplorer({ root, repository });
// An explicit instance also makes local reproducibility and accessibility checks possible.
globalThis.rnaExplorer = explorer;
explorer.start().catch(error => { explorer.status(`RNA data could not be loaded: ${error.message}`, 'error'); console.error(error); });
