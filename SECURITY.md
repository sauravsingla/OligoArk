# Security Policy

## Supported versions

OligoArk is pre-1.0 research software. Security fixes are applied to the latest `main` branch.

## Automated checks

Pull requests and `main` are checked with Ruff, strict mypy, Bandit, pip-audit, package installation checks, Docker build validation, and GitHub CodeQL. Dependabot monitors Python and GitHub Actions dependencies. These scanners are review aids, not evidence that the reference API is production hardened.

## Reporting

Please use GitHub private vulnerability reporting when it is enabled for the repository. Do not place sensitive production data, credentials, private sequencing data, or exploit details in a public issue.

## Current hardening boundaries

- `DNAArchive` validates archive format, metadata types, strand alphabet, declared counts, and checksum metadata, but the library intentionally does not impose a universal file-size limit.
- The reference FastAPI service enforces `OLIGOARK_MAX_API_PAYLOAD_BYTES` on decoded input and serialized archive requests; deployments should also enforce reverse-proxy/body limits.
- Optimizer search and graph/alignment reconstruction can be computationally expensive on adversarial or very large inputs. Production deployments need quotas, timeouts, and rate limiting.
- The reference API has no authentication or authorization and is intended for local/research deployment unless placed behind an authenticated gateway.
- The provided Docker image runs as an unprivileged user, but container isolation is not a substitute for application-level access controls.
