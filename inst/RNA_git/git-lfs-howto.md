# Git LFS tips

Rstudio server doesn't use git lfs properly.  Some third-party GUIs for git seem to support it, such as [smartgit](https://www.syntevo.com/smartgit/download/), or gitkracken. (I can't specifically endorse either though).

## Workarounds to using LFS while staging in rstudio.

0.  Open a terminal via ssh or the rstudio in the project directory and make sure git and git-lfs are loaded:

        module load git/2.30.1 git-lfs/2.13.2
        eval `ssh-agent`; ssh-add #optional, so you don't have retype your password constantly
1.  Stage / unstage files as normal in rstudio.
2.  Re-add the files (needed because lfs does something magical when git add is called). In your terminal:

        git add --renormalize .
3.  Check that git lfs is tracking properly.

        git lfs status
For the lfs files it should appear in this list and have `LFS: <hashcode>` next to the filename.  If not, make sure it's being tracked with `git lfs track <path/to/file.ext>`, and then rerun step 2.
4.  `git commit`
5.  `git push`

## Uhoh, something is broken:

1.  I see "File should be pointer" message when I pull.  This means that a file was added to the main repo rather than lfs.  We can fix this going forward by calling:

        git lfs migrate import --no-rewrite <pathspec>
