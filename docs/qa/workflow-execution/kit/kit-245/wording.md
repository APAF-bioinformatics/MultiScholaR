# KIT-245 Part S: drafts for the owner's wording (VAL-41)

Every user-facing text Part S adds, as written on `kit-245`. The owner sets the final words; nothing here depends
on OD-1 or OD-1b (the installed application's panel is unchanged).

## Messages (shown after "GitHub push failed: " in the Summary tab, or as the error in R)

| Class | Text | Where |
|---|---|---|
| `multischolar_github_project_id_invalid` | the Project ID may contain only letters, digits, '.', '-' and '_' (it becomes the repository's name on GitHub); nothing was created. | `R/func_github_management.R`, item 0 |
| `multischolar_github_git_missing` | pushing needs the git program, which was not found on this computer, so nothing was created on GitHub. Install Git (git-scm.com), start R again and press Push again. | item 2 |
| `multischolar_github_git_too_old` | pushing needs git 2.31 or later; this computer has &lt;version&gt;, so nothing was created on GitHub. Update Git (git-scm.com), start R again and press Push again. | item 1b |
| `multischolar_github_file_too_large` | &lt;file(s)&gt; is/are larger than GitHub's limit of 100 MB for one file, so nothing was created on GitHub. Remove it/them from &lt;scripts/omic&gt; or keep the project in Git another way. | item 3 |
| `multischolar_github_push_failed` | the repository &lt;org&gt;/&lt;id&gt; was created on GitHub but the files could not be sent (&lt;reason&gt;). It is empty: delete it on GitHub or use another Project ID. | item 4 |

## Help (`?pushProjectToGithub`, `?pushProjectToGithubFromDirs`)

Description of both: "Sends the omic's source folder as it is, including copies of input files and saved session
data."

Details of `pushProjectToGithub`: the paragraph in `man/pushProjectToGithub.Rd` (`\details`): the checks made before
anything is created on GitHub; the repository always private in the named organisation and what private means there;
what the commit holds per omic; the token never written to a file; the temporary copy removed whatever the outcome;
an empty repository named after a failure GitHub reports. `pushProjectToGithubFromDirs` points to it.

## NEWS (`### Fixes`, under the unreleased heading)

The KIT-245 bullet, and the qualification of KIT-234's 'Push to GitHub' bullet (P10: "It also needs the git program,
and it sends the omic's scripts folder, which holds copies of the input files (lipidomics, metabolomics,
transcriptomics) or the saved session with its data (proteomics)"), as in `NEWS.md` on `kit-245`.

## Product part (b) (VAL-244): the installed application's Summary tab

Under "Version Control (Optional)", in place of the checkbox, the four fields and the button
(`R/mod_omic_summary.R`, `mod_omic_summary_ui`):

> Pushing to GitHub is not available in the installed application. To keep this analysis under version control, use
> your usual Git tool; note that the project's scripts folder holds copies of your input files and analysis data,
> which you may not be allowed to share.

Difference from the ticket's draft (section 7): "the folder scripts/<omic>" is written "the project's scripts
folder", because the Summary tab's configuration does not carry the omic's folder name and adding it would change the
configuration schema (`R/utils_omic_gui_config.R:377`). The owner may prefer the folder named; that is a small
follow-up.

NEWS: the KIT-245 (b) bullet under Fixes, and KIT-234's 'Push to GitHub' bullet in P10's (b) form.

## Open for the owner

- Whether "private means visible to that organisation's owners and to the members its settings allow" is the right
  plain description (GitHub's model, not re-read; no external call).
- The R users' message with no token ("DEBUG: Failed to retrieve GitHub token from gh::gh_token()") is unchanged
  (out of scope, F6).
