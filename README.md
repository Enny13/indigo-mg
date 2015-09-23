indigo-mg - Indigo Open Molecule Generator

Introduction
========

You are currently reading the README file for the **indigo-mg** Project. 
This project is based on ideas of **OMG** project which is hosted on http://sourceforge.net/p/openmg 

**OMG** (Copyright 2011 The Netherland Metabolomics Center Team License: GPL v3, see License-gpl-3.txt) is an open-source tool for the generation of chemical structures, implemented in the programming language Java(tm).
The library is published under terms of the standard The GNU General Public License (GPL) v3.

**indigo-mg** replaces **CDK** library which serves as tool for molecules instrumentation in **OMG** with **Indigo Toolkit**.
Other major difference is that methods for generation of molecules lists are separated from writing to output files. It makes easier to use this code as a library.   
        
**indigo-mg** is at the moment about 2 times slower than **OMG**, but optimization-refactoring is in progress.

PLEASE NOTE: **indigo-mg** is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Important note
=======

Note that valence model in **CDK** and in **Indigo Toolkit** differs, so the count of generated molecules differ in **OMG** and **indigo-mg** 

System Requirements
========

Java 8 is required.
Maven 3.0.5 or higher 

Using indigo-mg tool
==============

First ```MoleculesGenerator``` should be created from ```Indigo``` instance:
```java
Indigo indigo = new Indigo();
MoleculesGenerator gen = new MoleculesGenerator(indigo);
```
After that you can use one of methods called ```generateMolecules``` which differ by format in which elemental composition and substructures are given. For example,
    
```java
List<IndigoObject> mols = gen.generateMolecules("C3H6", "C(O)=O");
```

Here elemental composition is given by string of format {elname_1}{elcount_1}{elname_2}{elcount_2}...{elname_n}{elcount_n} and substructures are given by smiles strings separated with whitespace.

For using tool as standalone application, see [Setup and Build section](#setup) 

File structure
==============

1. src/ : Contains source code of the application
2. test/ : Contains tests of application
3. pom.xml : Maven pom file
4. License-gpl-3.txt : The license file.
5. README.md : This file.

Setup and build<a name="setup"></a>
===============

To build executable jar from the root folder of project, command

```bash
mvn package -P exec
```

The jar will be at target folder. Now you can generate molecules with command

```bash
java -jar indigo-mg-1.0-jar-with-dependencies.jar -ec C6H8 -o out.sdf
```

In out.sdf you will have all molecules with given element composition. If you like to specify substructures (only non-overlapping are allowed now), you can provide sdf file with them, say subs.sdf: 
 
```bash
java -jar indigo-mg-1.0-jar-with-dependencies.jar -ec C6H8 -fr subs.sdf -o out.sdf
```




