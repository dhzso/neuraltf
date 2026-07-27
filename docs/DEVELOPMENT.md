# Development Guide

> Engineering standards and workflow for BioForge contributors.

---

## Repository Workflow

1. **Read the specification.** Every task begins with a written specification
   that defines scope, deliverables, and acceptance criteria.

2. **Implement only the requested scope.** Do not continue into future
   milestones. Do not add unrequested files or features.

3. **Explain what was changed.** Provide a summary of files created and
   modified.

4. **Explain why.** Provide rationale for implementation decisions.

5. **Report assumptions.** Document any assumptions made during implementation.

6. **Do not continue into future milestones.** Stop after completing the
   current task.

---

## Commit Message Convention

Follow the [Conventional Commits](https://www.conventionalcommits.org/) format:

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

### Types

| Type       | Purpose                                           |
|------------|---------------------------------------------------|
| `feat`     | A new feature                                     |
| `fix`      | A bug fix                                         |
| `docs`     | Documentation-only changes                        |
| `style`    | Code style changes (formatting, whitespace)       |
| `refactor` | Code restructuring without behavior change        |
| `test`     | Adding or modifying tests                         |
| `chore`    | Maintenance tasks (build, dependencies, config)   |
| `ci`       | CI/CD pipeline changes                             |

### Examples

```
docs: add repository foundation documentation
feat(scrna): add quality control filtering module
fix(docker): pin base image to debian-slim-12.5
```

---

## Branch Strategy

### Branch naming

```
<type>/<scope>-<short-description>
```

Examples:

```
docs/repository-foundation
feat/scrna-quality-control
fix/docker-base-image
```

### Branch model

- `master` — stable, always in a working state
- Feature branches — created from `master`, merged back via pull request
- No direct commits to `master` except for the initial repository setup

---

## Coding Standards

| Standard           | Requirement                              |
|--------------------|------------------------------------------|
| Language           | Python 3.11                              |
| Style              | PEP8                                     |
| Type hints         | Required on all function signatures      |
| Comments           | Meaningful and only where necessary      |
| Imports            | Explicit, no wildcard imports            |
| Line length        | 88 characters (Black default)            |
| Naming             | `snake_case` for functions and variables |
| Classes            | `PascalCase`                             |
| Constants          | `UPPER_SNAKE_CASE`                       |

### Additional rules

- Never use "latest" image tags — always pin explicit versions
- Never introduce unnecessary dependencies
- Never delete datasets or move folders without approval
- Never silently rename files

---

## Versioning

BioForge follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html):

```
MAJOR.MINOR.PATCH
```

During the alpha phase, versions use the pre-release suffix:

```
0.1.0-alpha.1
```

The changelog follows the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format.

---

## Code Review Process

1. **Open a pull request** from the feature branch to `master`.

2. **Self-review** the diff before requesting review. Verify:
   - Only requested scope is implemented
   - No unrequested files are added
   - Code meets the coding standards above
   - Commit messages follow the convention

3. **Request review** from the Architecture Lead or Owner.

4. **Address feedback** with additional commits on the same branch.

5. **Merge** only after approval. Squash-merge is preferred to keep the
   commit history clean.

6. **Delete the feature branch** after merge.
