# Contributing to 3D-POLE-Render

Thank you for your interest in contributing! This guide explains how to get started.

## Getting Started

1. Fork the repository and clone your fork.
2. Install dependencies:
   ```bash
   npm install
   ```
3. Run the linter and tests to confirm everything works:
   ```bash
   npm run lint
   npm test
   ```

## Development Workflow

1. Create a feature branch from `main`:
   ```bash
   git checkout -b your-feature-name
   ```
2. Make your changes, following the coding standards below.
3. Add or update tests for any new functionality.
4. Ensure linting and tests pass:
   ```bash
   npm run lint
   npm test
   ```
5. Commit your changes with a DCO sign-off (see below).
6. Open a pull request against `main`.

## Coding Standards

- **JavaScript**: Follow the ESLint configuration in `eslint.config.js`. Run `npm run lint:fix` to auto-fix issues.
- **Tests**: Write tests using [Vitest](https://vitest.dev/). Place test files in `test/js/`.
- **Julia pipeline**: Follow existing conventions in `pipeline/src/`. Run the Julia test suite with `julia --project=pipeline -e 'using Pkg; Pkg.test()'`.

## Pull Request Guidelines

- Keep PRs focused — one feature or fix per PR.
- Write a clear description of what changed and why.
- Reference any related issues (e.g., `Fixes #42`).
- Significant changes (new features, architecture, dependencies) require approval from the Project Lead per the [Governance](GOVERNANCE.md) policy.

## Developer Certificate of Origin (DCO)

All contributions must include a DCO sign-off. This certifies that you have the right to submit the contribution under the project's license. Add the following to each commit message:

```
Signed-off-by: Your Name <your.email@example.com>
```

You can do this automatically with:

```bash
git commit -s -m "Your commit message"
```

By signing off, you agree to the [Developer Certificate of Origin](https://developercertificate.org/).

## Reporting Issues

- Use [GitHub Issues](https://github.com/bloomed-health/3D-POLE-Render/issues) for bugs and feature requests.
- For security vulnerabilities, see [SECURITY.md](SECURITY.md).

## License

By contributing, you agree that your contributions will be licensed under the [Apache 2.0 License](LICENSE).
