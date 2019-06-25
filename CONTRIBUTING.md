# DEVELOPER GUIDE and WISH LIST
## CONTRIBUTING
You are very welcome to give your contribution to **thermocepstrum**.

Here are a few guidelines to consider when contributing.

- The **`master`** branch is protected and contains the latest stable releases only. Releases are tagged, e.g. `v0.1.2`.
- The **`develop`** branch hosts the main development branch. Pull requests from forks/branches should be made to this branch.
- If you start working on a **new feature**, you should open a *new* branch (from `develop`), labelled with the name of the feauture you want to implement.
- If you find any **issue** or want to suggest a new feature/fix that can be implemented in the future, open a new *Issue* describing the problem. In this way other users can help or give their contribution.
- If you want to implement a **solution** to an issue, name that branch as `fix_ISSUENUMBER_ISSUENAME`.
- When you are done with a branch and want to **merge** it to `develop`, submit a *Pull request*. Unless it is a simple fix, and especially if it is a new functionality, this will help to keep the code clean and find potential bugs beforehand. Also, it will help to standardize the style.
- Before updating the `develop` branch, make **tests**!
- Try to make **commits** as specific as possible, and document the changes.


## WISH LIST

### WISH LIST FOR v0.2/3 RELEASE:
- [x]  write standard setup.py and package info files
- [ ]  i/o example
- [x]  move `grafici_belli` to utils


### WISH LIST FOR v1.0 RELEASE:
#### easier stuff:
- [ ] Python 3 migration
- [ ] Nice formatting
- [ ] translate all comments to English
- [ ] renew plot style
- [ ] i/o: when reading a text file, the program should save a binary file (unless otherwise stated). If it finds the binary file it does not read the text one (making it clear to the user)
- [ ] in plot functions, always return (figure, axes)

#### harder stuff:
- [ ] Rethink about base data classes: MDSample should be the most general, HeatCurrent should be a subclass specific for thermal transport. In future we may want to extend it to other transport coefficients.
- [ ] CepstralResults class (to discuss)
- [ ] Cepstral class (subclass of CosFilter)
- [ ] Plot class that takes care of the plots from a CepstralResults object or MDSample (to discuss)
- [ ] thermocepstrum GUI
- [ ] there should be one single `compute_periodogram` method (either for single and multi component)

