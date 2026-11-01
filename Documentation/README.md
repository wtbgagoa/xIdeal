# xIdeal documentation

The editable Wolfram Documentation Tools source notebooks are in:

- `English/ReferencePages/Symbols/`: one page for each of the symbols with an explicit usage message in the kernel source.

The pages include usage syntax, implementation notes, See Also links, package initialization, self-contained example setup, evaluated examples, and categorization metadata. The documentation link base is `xIdeal`.

## Editing

Open a source notebook in Mathematica and use **Palettes > Documentation Tools**. Edit these source notebooks, not the generated notebooks under `build/`. Run each page's example initialization and setup once in a fresh kernel with xAct installed. 

## Building and previewing

Use a Mathematica installation with PacletTools, Documentation Tools, and a working notebook front end. xAct is needed to evaluate the examples. In a fresh kernel, set `pacletRoot` to the directory containing `PacletInfo.wl`, then evaluate:

```wolfram
Needs["PacletTools`"];
result = PacletBuild[pacletRoot];
```

`PacletBuild` compiles the documentation and creates `build/xTerior-0.1.0.paclet`. For a documentation-only build, use `PacletDocumentationBuild[pacletRoot]`. Compilation is necessary: simply archiving the editable source notebooks does not produce finished Documentation Center pages.

Preview the built pages without replacing an installed paclet:

```wolfram
PacletDirectoryLoad[FileNameJoin[{pacletRoot, "build", "xIdeal"}]];
SystemOpen["paclet:xTerior/ref/xActSolutionData"];
```

To install the built archive, evaluate:

```wolfram
PacletInstall[result["PacletArchive"], ForceVersionInstall -> True]
```

`Get[FileNameJoin[{pacletRoot, "makexIdeal.wl"}]]` performs the build and installation together. A failed build leaves the installed paclet untouched. The package version remains 0.1.0; increment it before distributing a new release if appropriate.

The file layout, source stylesheets, metadata, and build process follow Wolfram's [Creating Paclets](https://reference.wolfram.com/language/PacletTools/tutorial/CreatingPaclets.html), [Authoring Tech Notes](https://reference.wolfram.com/language/DocumentationTools/tutorial/AuthoringTechNotesUsingDocumentationTools.html), and [PacletDocumentationBuild](https://reference.wolfram.com/language/PacletTools/ref/PacletDocumentationBuild.html) documentation.
