# CLAUDE.md guidance for epib705 repo

---

## Project Overview

- This repo contains slides for delivering a lecture to PhD students in epidemiology on the concepts and methods behind quantitative bias analysis. 
- It draws heavily on the book by Fox, MacLehose, and Lash (2021) Applying Quantitative Bias Analysis to Epidemiologic Data. You can find a copy of the book at "/Users/samharper/Library/CloudStorage/Dropbox/work/research/zotero/Fox-2021-Applying.pdf"

- The purpose is to show students the basics of why you may want to do quantitative bias analysis (QBA), some simple examples of how to do it, and to give them some additional resources.
- The lecture is designed for a two-hour class, with a break in the middle. 
- The basic format is a quarto presentation (.qmd) that will be rendered to .html for in-class presentation but can later be exported to .pdf, so the slides should look good for both formats.

---

## Structure of the Repo

- Keep all images in `images/`. Any images you may create or download from elsewhere should be store here. 

- If there is any analysis (tabular, regression, etc.) of .csv or other data files, those files should be store in a folder called `data/`. 
- References in the slides are given at the end and are contained in the `bias-analysis.bib` file. You can add to this file as needed, but make sure to use the same format for consistency.

- The basic format for the slides is in `custom.scss`.

## Instructions for Claude

- I want you to update the slide deck for 2026 and create a new .qmd file (the last version is shown as `EPIB705-Bias-Analysis-2024.qmd`). The new version should be called `EPIB705-Bias-Analysis-2026.qmd`. DO NOT overwrite the existing 2024 version, as I want to keep that for reference.
- Read through the existing .qmd file as well as the rendered .html version to get a sense of my slide style, proportions, white space, etc.
- I want you to try and clean up the .qmd file and make it as readable as possible without sacrificing content or completeness. Try and reduce unnecessary spaces and tabs. Make the visuals as clean and simple as possible. Try and center tables where possible across the whole slide or within a column. 
- One of the main challenges with the existing deck is that to include variables and equations in the data tables I current have to use native latex code, like this: 

```{=latex}
\begin{tabular}{cccccccccccc}
\hline 
 & \multicolumn{3}{c}{{\small{}Z=1 Smokers}} &  & \multicolumn{3}{c}{{\small{}Z=0 Non-smokers}} &  & \multicolumn{3}{c}{{\small{}Total}}\tabularnewline
 & {\small{}X=1} & {\small{}X=0} & {\small{}Total} &  & {\small{}X=1} & {\small{}X=0} & {\small{}Total} &  & {\small{}X=1} & {\small{}X=0} & {\small{}Total}\tabularnewline
\hline 
{\small{}D=1} & {\small{}${\scriptstyle A_{11}}$} & {\small{}${\scriptstyle A_{01}}$} & {\small{}${\scriptstyle M_{11}}$} &  & {\small{}${\scriptstyle A_{1+}-A_{11}}$} & {\small{}${\scriptstyle A_{0+}-A_{01}}$} & {\small{}${\scriptstyle M_{A+}-M_{11}}$} &  & {\small{}${\scriptstyle A_{1+}}$} & {\small{}${\scriptstyle A_{0+}}$} & {\small{}${\scriptstyle M_{A+}}$}\tabularnewline
{\small{}D=0} & {\small{}${\scriptstyle B_{11}}$} & {\small{}${\scriptstyle B_{01}}$} & {\small{}${\scriptstyle M_{01}}$} &  & {\small{}${\scriptstyle B_{1+}-B_{11}}$} & {\small{}${\scriptstyle B_{0+}-B_{01}}$} & {\small{}${\scriptstyle M_{B+}-M_{01}}$} &  & {\small{}${\scriptstyle B_{1+}}$} & {\small{}${\scriptstyle B_{0+}}$} & {\small{}${\scriptstyle M_{B+}}$}\tabularnewline
\hline 
{\small{}Total} & {\small{}${\scriptstyle N_{11}}$} & {\small{}${\scriptstyle N_{01}}$} & {\small{}${\scriptstyle N_{+1}}$} &  & {\small{}${\scriptstyle N_{1+}-N_{11}}$} & {\small{}${\scriptstyle N_{0+}-N_{01}}$} & {\small{}${\scriptstyle N_{++}-N_{+1}}$} &  & {\small{}${\scriptstyle N_{1+}}$} & {\small{}${\scriptstyle N_{0+}}$} & {\small{}${\scriptstyle N_{++}}$}\tabularnewline
& & & & & & & & & & & \tabularnewline 
\end{tabular}
```

- This is a huge pain and I hope there is some way to make this easier because it is really a pain to write all of that code to make a table and to update it.

---

## Table Sizing and Font Control

- When a `tinytable` (`tt()`) table is too wide or tall to fit comfortably on a slide, **do not** add `{.smaller}` to the slide heading. That class shrinks *all* text on the slide — headings, bullet points, and equations — making the presentation harder to read.
- Instead, shrink only the table font by passing `fontsize` to `style_tt()`. A value of `0.9` (i.e., 90% of the default font size) is usually sufficient for large tables:

```r
tt(my_data, escape = FALSE) |>
  group_tt(j = list("Group A" = 2:4, "Group B" = 6:8)) |>
  style_tt(j = 1:8, fontsize = 0.9)
```

- Apply `fontsize` to all columns (`j = 1:ncol`) so the header row is also scaled consistently.
- `fontsize` values below `0.7` are generally too small to read comfortably in a presentation; prefer restructuring the table or splitting it across slides instead.
- The `{.smaller}` slide class should be avoided entirely in this deck.

---

## TikZ / DAG Figure Sizing

- The YAML header sets `auto-stretch: false` under `format: revealjs:`. This disables Quarto's default behaviour of automatically applying the `r-stretch` CSS class to figures, which would otherwise make a DAG expand or shrink to fill whatever vertical space remains on a slide (causing the same diagram to appear at wildly different sizes depending on how much other content is present).
- Control TikZ figure size using `fig-height` (a number in inches) directly in the chunk options. Do **not** use `fig-width: <percentage>` — that is invalid YAML and will cause a render error. Use `out-width: <percentage>` if you need to constrain display width, but prefer `fig-height` for TikZ.
- For the simple three-node horizontal DAGs in this deck (e.g. Z → X → D), `fig-height: 1.5` produces a compact, consistently-sized figure that leaves room for tables and bullets below it. Adjust up or down to taste, but keep the value the same across all equivalent DAG slides so they look uniform.

```r
# Correct — fig-height is a number (inches), consistent across DAG slides
```{r, engine='tikz'}
#| fig-align: center
#| fig-height: 1.5
\begin{tikzpicture}...
\end{tikzpicture}
` ` `
```
