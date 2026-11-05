# Nuova architettura degli stati

Questa directory contiene l’implementazione parallela della gestione degli
stati. Non viene ancora usata dai consumatori di produzione e non deve
dipendere da `state_conversion_2d`.

Il nucleo è organizzato attorno a prodotti computazionali distinti:

- `conservative_to_primitive`
- `primitive_to_conservative`
- `evaluate_mixture`
- `build_physical_state_from_q`
- `build_physical_state_from_qp`
- `update_physical_state_kinematics`

Le formule REAL/COMPLEX comuni sono mantenute in include `.inc`; i wrapper
tipati forniscono soltanto dichiarazioni, workspace e gestione dei branch
necessari al complex-step. Sono ora presenti i kernel di valutazione della
miscela e le conversioni conservative/primitive, ancora non collegati ai
consumatori di produzione.

La migrazione dei consumatori e l’inserimento nel target principale inizieranno
solo dopo i test di caratterizzazione e il confronto con l’implementazione
attuale.
