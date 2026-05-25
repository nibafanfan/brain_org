# LLM prompt — expand collapsed Samples column into per-GSM annotation sheet

Paste the contents below into a fresh LLM session, attach `brain_organoid (3).xlsx`, and let it produce a new annotation workbook.

---

## Role

You are a **bioinformatics data scientist** working on a brain-organoid scRNA-seq atlas. Your task is to convert a paper-level annotation spreadsheet (one row per GEO deposit) into a per-GSM annotation spreadsheet (one row per individual sample), matching the structure of an existing reference workbook.

## Input

`brain_organoid (3).xlsx`, sheet **"Paper Datasets"**. Each data row (row 2 onward) is one GEO deposit, with these columns:

| Col | Name | Notes |
|---|---|---|
| 1 | Flag | Usually empty |
| 2 | Search Modality | How the deposit was discovered (`Self`, `GPT`, `GEO`) |
| 3 | Paper | Publication title |
| 4 | Dataset | A hyperlink cell — the URL behind it points to the GEO accession (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE######`). **Read the hyperlink target, not the cell text** — the cell often renders as `"GEO Accession viewer"` |
| 5 | (unnamed extra hyperlink) | Sometimes used; same treatment |
| 6 | Dataset Notes | e.g., `FASTQ Only`, `NOT PRIMARY DATA` |
| 7 | Cell Type Origin | e.g., `"HMGU1; Human induced pluripotent stem cells"`, `"H9; Human embryonic stem cells"` |
| 8 | Age/Timepoint | e.g., `"69 Days"`, `"30 days; 90 days; 150 days"`, `"1 month"` |
| 9 | Guided vs Unguided | `Guided` or `Unguided` (sometimes multiline) |
| 10 | Organoid Type | e.g., `Cerebral`, `Midbrain`, `Cortical, Thalamic, Spinal, Sensory` |
| 11 | Protocol | e.g., `Conventional`, `Vascularized by transgene`, `Conventional + Microglia`, `Air-Liquid Interface Slice Organoids` |
| 12 | Multi-Lineage | `Yes` or `No` |
| 13 | **Samples** | Multi-line text — see filter and parse rules below |
| 14 | Notes | Free-form |

## Filter rule (which deposits to include)

Process **only** the deposits whose `Samples` column (col 13) is **one of**:

- **(A) Explicit GSM list** — cell contains one or more `GSM\d+` IDs (each on its own line, usually accompanied by a sample label).
- **(B) Exactly "all samples"** — cell text, stripped and lowercased, equals exactly `"all samples"` (no qualifying suffix, no extra lines).

**Skip** all other rows:
- `"All samples - xxx"` or `"All samples (VEH ONLY)"` or `"All samples\nWT and Mutant are mixed..."` etc. — these are qualified and need colleague disambiguation.
- Paper-specific descriptive text without GSM IDs (e.g., `"WT Samples Only"`, `"Control (CS01ictr)"`).
- Empty Samples cells.
- Rows whose Dataset column has no GEO hyperlink (no parseable accession).

## Output

A single Excel workbook `brain_organoid_annotation_expanded.xlsx` with one sheet `"GSM Annotations"`, one row per GSM, matching the schema of `brain_organoid_annotation.xlsx`'s `"GSM Annotations"` sheet:

| Col | Name | How to fill it |
|---|---|---|
| 1 | Paper | Copy from input col 3 |
| 2 | GEO | The GEO accession (`GSE######`) extracted from the Dataset hyperlink target |
| 3 | GSM | The individual GSM ID |
| 4 | GSM Label | The descriptive text on the same Samples-cell line as the GSM, with the `GSM\d+` and any leading punctuation/whitespace stripped. Example: from `"GSM8501109 - WNTi (batch A)"` → `"WNTi (batch A)"`. For "all samples" rows, GSM Label is `"n/a"` (or fetched from GEO `Sample_title` if the LLM has GEO access). |
| 5 | Barcode ID | Default `"n/a"`. (Filled only when the colleague has noted a per-sample barcode prefix — leave `n/a` unless input gives one.) |
| 6 | Cell Type Origin | Normalize the input col 7. Rule: if it contains `"induced pluripotent"` → `"ipsc"`. If it contains `"embryonic"` → `"esc"`. Drop the cell-line prefix (e.g., `"HMGU1; Human induced pluripotent stem cells"` → `"ipsc"`). For mixed-origin deposits, use the line corresponding to this sample if disambiguable, else use the first-listed origin. |
| 7 | Age (Days in Vitro) | Parse input col 8 to an integer day count. Conversions: `N days` → `N`, `N weeks` → `N*7`, `N months` → `N*30`, `D60` → `60`. If the cell has multiple ages (one per timepoint, joined by `\n` or `;`), match the age to this specific GSM if disambiguable from the GSM Label — otherwise leave as the first age value. |
| 8 | Organoid Type | Copy from input col 10. If multi-region (`"Cortical, Thalamic, Spinal, Sensory"`), keep as comma-separated. |
| 9 | Protocol | Copy from input col 11. |
| 10 | Unguided? | `Yes` if input col 9 contains `Unguided`, else `No`. |
| 11 | Multi-Lineage? | Copy `Yes`/`No` from input col 12. |
| 12 | Vascularized? | `Yes` if input col 10 (Organoid Type), col 11 (Protocol), or col 13 (Samples label for this GSM) mentions any of: `Vascular`, `vasculariz`, `endothelial`, `pericyte`, `perfusable`. Else `No`. |
| 13 | Slice? | `Yes` if input col 11 (Protocol) or the GSM Label contains `slice` (case-insensitive, whole-word) or `air-liquid` or `ALI`. Else `No`. |

## Expansion rule for "all samples" deposits

When the Samples cell is exactly `"all samples"`:

1. Use the GEO Soft-format API to enumerate every GSM under that accession:
   - URL: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>&targ=self&form=text&view=brief`
   - Parse lines starting with `!Series_sample_id = ` to get the GSM list.
2. For each GSM, fetch its individual `Sample_title` via:
   - URL: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSM>&targ=self&form=text&view=brief`
   - Line: `!Sample_title = "<title>"`
3. Emit one output row per GSM with `GSM Label` = the fetched Sample_title.

If GEO access is unavailable, leave `GSM Label = "n/a"` and skip the per-GSM age disambiguation step (use the deposit-level age value).

## Edge cases to handle

- **Multiline Samples cells with mixed control/treated GSMs already filtered**: trust the xlsx — anything in the cell is a control GSM by design.
- **Tab vs space between GSM and label**: both are valid separators in the Samples cell. Treat any whitespace as the separator.
- **Duplicate GSMs across multiple lines within one cell**: legitimate if the surrounding labels differ (e.g., `C1` vs `C2` conditions). Emit one row per *unique (GSM, label)* pair; if both occurrences have identical labels (whitespace-only difference), emit just one row.
- **Reference rows** (Search Modality = `Self`, no GEO link): if Dataset col 4 has no GSE accession, skip the row.
- **Multiple GSE in one hyperlink**: if the Dataset URL contains multiple accessions (rare), use only the first; flag in a `Notes` column.

## Quality checks before writing the output

1. Row count = total individual GSMs emitted (expect ~400–700 depending on how many "all samples" deposits expand).
2. No duplicate (GEO, GSM, GSM Label) triples.
3. Every `Age (Days in Vitro)` is a positive integer or empty.
4. Every `Unguided?`, `Multi-Lineage?`, `Vascularized?`, `Slice?` is `Yes` or `No`.
5. Every `Cell Type Origin` is `ipsc`, `esc`, or `n/a`.

## Output file

`brain_organoid_annotation_expanded.xlsx` with sheet `"GSM Annotations"` matching the 13-column schema above. Include a brief stats summary as a second sheet `"Summary"`:

- Total deposits processed
- Total deposits included via "all samples" path
- Total deposits included via explicit-GSM path
- Total GSMs emitted
- Breakdown of `Vascularized? = Yes`, `Multi-Lineage? = Yes`, `Slice? = Yes`, `Unguided? = Yes`
- Median `Age (Days in Vitro)` and range

## Style

- Use openpyxl or pandas. Don't reformat colors or sheet layout — match the reference `brain_organoid_annotation.xlsx` cell-style for the `GSM Annotations` sheet.
- Don't add columns beyond the 13 defined.
- Don't filter or rank rows beyond what this prompt specifies.
- Preserve the original Paper title text verbatim (don't truncate).
