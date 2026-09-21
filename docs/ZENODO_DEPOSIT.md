# Minting a DOI for this archive

Reference 11 of the manuscript cites this repository. Every reviewer so far has
said a version tag alone does not freeze the submission if the files cannot be
retrieved, and the repository is currently **private**. This is the step that
unblocks that.

Two routes. Route A is better if the repository is going public anyway. Route B
mints a DOI while the GitHub repository stays private.

---

## Before either route: three things to decide

**1. The licence you select on Zenodo is not the licence of this archive.**
Zenodo asks for one licence per record. This archive has two:

- analysis code: MIT (see `LICENSE`)
- AlphaGenome predictions and everything derived from them: DeepMind's
  AlphaGenome Output Terms of Use, which are **non-commercial** (see
  `LEGALLY_BINDING_TERMS_OF_USE.txt`)

Selecting a blanket MIT or CC-BY on the Zenodo record would assert rights over
the model outputs that we do not hold. **Choose a non-commercial or "other"
licence** and state the split in the description. Suggested description text is
at the bottom of this file.

**2. The record must be Open Access**, not restricted. A restricted record
reproduces exactly the problem the reviewers raised.

**3. Confirm nothing sensitive is being published.** Already checked on
21 September 2026: no API key in the working tree or in git history, and
`analysis/ag_auth.py` reads the key from `~/.alphagenome/api_key`, outside the
repository. Re-run before publishing:

    grep -rInE "AIza[0-9A-Za-z_-]{20,}|sk-[A-Za-z0-9]{20,}" --exclude-dir=.git .
    git log --all -p | grep -InE "AIza[0-9A-Za-z_-]{20,}|sk-[A-Za-z0-9]{20,}"

Both should return nothing.

---

## Route A: GitHub release, archived automatically

**The order matters. Zenodo only archives releases created AFTER the repository
is switched on. A release made first is not picked up.**

1. Make the repository public:
   `gh repo edit PatientThread/apol1-alphagenome --visibility public`
   (or Settings, Danger Zone, Change visibility)
2. Sign in at https://zenodo.org with the GitHub account that owns the repo.
3. Go to https://zenodo.org/account/settings/github/ and toggle
   **PatientThread/apol1-alphagenome** to ON. If it is not listed, use Sync now.
4. Create a **release**, not just a tag. The tag `v1.0-bmcrn` already exists:

       gh release create v1.0-bmcrn \
         --title "v1.0-bmcrn: BMC Research Notes submission archive" \
         --notes "Frozen analysis code and results supporting the manuscript."

5. Zenodo receives the webhook within a few minutes and mints **two** DOIs:
   - a **concept DOI**, which always resolves to the newest version
   - a **version DOI**, which resolves to this release only
   **Cite the version DOI in the manuscript.** The concept DOI would let the
   cited content change after review.
6. Open the new record, click Edit, and correct the metadata Zenodo guesses:
   author name and ORCID (0000-0002-8159-0879), title, licence (see above),
   and the description. Save and publish.

## Route B: manual upload, repository stays private

1. Build a clean archive from the tagged commit, so nothing untracked leaks in:

       cd "<repo>"
       git archive --format=zip --prefix=apol1-alphagenome-v1.0-bmcrn/ \
         -o ~/Desktop/apol1-alphagenome-v1.0-bmcrn.zip v1.0-bmcrn

2. Go to https://zenodo.org/uploads/new and upload that zip.
3. Set: Upload type **Software** (or Dataset), title, author with ORCID,
   description, keywords, version `v1.0-bmcrn`, and the licence discussed above.
4. Under Related identifiers, add the GitHub URL as "is supplement to" if you
   wish, though a private URL will not resolve for readers.
5. Publish. One DOI is issued; Zenodo also creates a concept DOI for future
   versions.

Note: the archive is about 390 MB, most of it ENCODE peak files. That is well
inside Zenodo's 50 GB limit. If you would rather keep it small, the ENCODE beds
can be omitted because `results/peak_manifest.csv` lists every accession and
the scripts re-download them, but keeping them makes the analysis reproducible
offline and I would keep them.

---

## After the DOI exists

1. Put the **version** DOI into reference 11 of the manuscript, replacing the
   GitHub release URL, and into `CITATION.cff` as `doi:`.
2. Rebuild the submission file:

       cd publication
       python3 <path>/md_to_docx.py MANUSCRIPT.md
       python3 format_for_submission.py

3. Check the DOI resolves in a private browser window before submitting.

## Licence text to paste

Zenodo's "Other" licence option takes a free-text description. The exact text
is in `docs/ZENODO_LICENCE_TEXT.txt` as a single unwrapped paragraph.

**Copy it from that file, not from this one.** An earlier version of this guide
presented the text as a markdown blockquote, and the leading "> " characters
were pasted into the live record along with it, where they display literally.
