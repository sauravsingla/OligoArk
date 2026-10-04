# Security Policy

## Supported versions

The project is pre-1.0 research software. Security fixes are applied to the latest `main` branch.

## Reporting

Please use GitHub private vulnerability reporting when enabled. Do not include sensitive production data in a public issue.

## Current hardening boundaries

- Archive JSON is size-unbounded at the library layer; services should enforce request/body limits.
- DNA decoding and edit-distance work can be computationally expensive on adversarial inputs.
- The reference API has no authentication or authorization and is intended for local/research deployment unless wrapped by a hardened gateway.
