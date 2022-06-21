
# Contributing to SPM.jl

Welcome to the `SPM.jl` repository, and thank you for thinking about contributing!

## Be sure to check out our Community Guidelines

The `SPM.jl` community follows 
[this Code of Conduct](https://github.com/spm/SPM.jl/blob/main/.github/CODE_OF_CONDUCT.md)
to ensure that our online spaces are enjoyable, inclusive, and productive for
all contributors.

## Developer guidelines

These guidelines are designed to make it as easy as possible to get involved.
If you have any questions that are not discussed in our documentation, or it is
difficult to find what you are looking for, please let us know by opening an
[issue](https://github.com/spm/SPM.jl/issues).

Please have a look at [Julia's Contributing Guidelines](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md).

### Coding Style

We will describe here the recommended coding style for the `SPM.jl` repository but in the meantime have a look at 
[Julia's Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/) and
[Invenia's Blue Style](https://github.com/invenia/BlueStyle). See also [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl).

In short the general formatting guidelines are:

* Use 4 spaces per indentation level, no tabs.
* Use whitespace to make the code more readable.
* No whitespace at the end of a line (trailing whitespace).
* Comments are good, especially when they explain the algorithm.
* Try to adhere to a 92 character line length limit.
* Use upper camel case convention for modules, type names.
* Use lower case with underscores for method names.
* Import modules with using, with one module per line and at the top of the file when possible.
* Avoid padding brackets with spaces. ex. `Int64(value)` preferred over `Int64( value )`.

### Pull Requests

Git recommendations for pull requests:

* Avoid working from the `main` branch of your fork, creating a new branch will make it easier if `SPM.jl`'s `main` changes and you need to update your pull request.
* Try to squash together small commits that make repeated changes to the same section of code so your pull request is easier to review, and `SPM.jl`'s history won't have any broken intermediate commits. A reasonable number of separate well-factored commits is fine, especially for larger changes.
* If any conflicts arise due to changes in `SPM.jl`'s `main`, prefer updating your pull request branch with `git rebase` versus `git merge` or `git pull`, since the latter will introduce merge commits that clutter the git history with noise that makes your changes more difficult to review.
* Descriptive commit messages are good.
* Using `git add -p` or `git add -i` can be useful to avoid accidentally committing unrelated changes.
* GitHub does not send notifications when you push a new commit to a pull request, so please add a comment to the pull request thread to let reviewers know when you've made changes.

#### Merges

Git recommendations for pull request reviewers:

* When merging, we generally like `squash+merge`. Unless it is the rare case of a PR with carefully staged individual commits that you want in the history separately, in which case `merge` is acceptable, but usually prefer `squash+merge`.

<!--

### Supported Versions

### Naming Conventions

### Testing

### Documentation

### Git Repository

#### Layout

#### Commit Messages

### Changelog

### Deprecations

-->