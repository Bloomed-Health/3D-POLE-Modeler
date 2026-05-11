# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability in this project, please report it responsibly. **Do not open a public issue.**

### How to Report

1. Email **security@bloomedhealth.com** with a description of the vulnerability.
2. Include steps to reproduce, affected versions, and any potential impact.
3. If possible, suggest a fix or mitigation.

### What to Expect

- **Acknowledgement** within 48 hours of your report.
- **Assessment** within 7 days — we will confirm the vulnerability and determine its severity.
- **Fix timeline** depends on severity:
  - **Critical/High**: Patch within 14 days.
  - **Medium/Low**: Patch in the next scheduled release.
- **Disclosure**: We will coordinate disclosure with you. We aim to publish a fix before any public disclosure.

### Scope

This policy covers the 3D-POLE-Render repository, including:

- JavaScript frontend code
- Julia pipeline code
- Build and CI configuration
- Documentation that could lead to misconfiguration

### Out of Scope

- Third-party dependencies (report upstream; we will update promptly)
- The RCSB PDB data service
- Issues in forked or modified versions of this project

## Security Best Practices for Contributors

- Never commit secrets, API keys, or credentials.
- Use HTTPS for all external data fetches.
- Keep dependencies up to date and review `npm audit` output.
- Follow the linting rules configured in `eslint.config.js`.

## Credits

We gratefully acknowledge security researchers who report vulnerabilities responsibly. With your permission, we will credit you in the release notes.
