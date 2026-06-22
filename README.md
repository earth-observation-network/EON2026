# EON Summer School website template

## Using the template

To use this template, you can click the "Use this template" button on the GitHub repository page, which will create a new repository with the same structure. 
**Remember to check the box to include all branches**.

## Updating the website

Please do not push new/updated .`qmd` files that run any code directly, but rather render the website locally (on your own computer), and then push the new files to GitHub.
To render the website locally, you need to use a terminal (not an R console):

```
quarto render
```

This will add new files (pre-rendered versions of the `.qmd` documents) to the `_freeze` folder.
You can then push these files to GitHub.
This push will automatically trigger a new GitHub Actions build, which will update the website.
