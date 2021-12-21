# {{Genesee.title}} {{Genesee.type}}

# Developer instructions

0. Consider having ssh-agent cache credentials.
    eval `ssh-agent`; ssh-add
1. Clone this repository
2. On bluehive, start rstudio-server
    rstudio-server.sh -S init-Rstudio-server.sh -- -t 0-6 --mem-per-cpu=32gb # modify node needs
3. proto_driverR.R is just a list of commands used to generate this repo
4. Modify scripts, copying changes that are generalizable back up to `{{Genesee.template.package}}/inst/{{Genesee.type}}_markdown` or add functions in `{{Genesee.template.package}}/R`
5. Big intermediate files that shouldn't be version-controlled can go in scratch/.  If they might be used by investigators, they should go in refined/ (which uses git LFS to avoid bloating the repo).
6.  Follow [Git-lfs instructions](git-lfs-howto.md) to commit into refined/ other large objects.  (It doesn't work from RStudio server, anyways.)


# Investigator instructions (WIP)

It is intended that this repo will serve as portal for investigators.

For single cell data, we'll describe

*  how to get the cloupe data
*  how to run iSEE locally
