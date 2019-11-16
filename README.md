# M-SPLIT
Notes on building this:

Using Eclipse for Java

1st step: Download jre 1.7/7 and then make eclipse workspace and import this repo from github

(import -> projects from git -> clone uri -> URL for this repo -> clone submodules checked)

right click msplit -> configure build path-> remove everything and add libraries->jre system library -> java se 7 (1.7)

then add dependencies (dependencies directory in this repo):

java build path -> libraries -> add external jar (for each jar in dependencies directory) and add external class folder for org directory in dependencies directory

change the library for all the projs to java 1.7 and fix the build paths (project -> properties -> java build path -> libraries tab->edit)

Make sure project build (project->properties->java build path->libraries tab-> JRE system library), compile (project->properties->java compiler->compiler compilance level), and run ((project->properties-> run/debug settings->edit->JRE tab) are all set to java 1.7 for EVERY project (or at least for the build and compile- change run for projects that we're running)

run MSPLIT/SpectrumLib to make SpectrumLib available under Launch configuration or MS/DecoySpectrumGenerator to make DecoySpectrumGene available for creating jar. to create jar, MSPLIT-> Export->Java-> Runnable JAR file -> extract required libraries
