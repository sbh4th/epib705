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
