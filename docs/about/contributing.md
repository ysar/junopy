# Contributing

JunoPy is an open-source package and anyone is welcome to create pull requests 
to add improvements. In order to prevent code bloat and retain ease of 
maintainability, only changes of the following types will be considered. 

1. *Fast*: Improvements to make code faster with minimal impact to readability 
or use.

2. *General*: Features that are general enough to be of use to a vast majority
 of users. If some analysis or processing task is specific to your research 
 and won't be used by others, please do not add that in.

3. *Simple*: Refactoring changes to reduce unnecessary or obsolete code. 

Please follow these guidelines:  

- Variables should be in `lowercase` or `camelCase`.  
- Class objects should be in `PascalCase`.  
- Methods/functions should be in `snake_case` with the appropriate verb as the
first word. E.g., `do_this_task`, `fetch_this_thing`, `calc_quantity`, etc.  
- Logicals should follow the variable naming style but with possessive verbs 
that describe the question. E.g. `hasData` or `isEmpty`. 
- Use descriptive but short names. If your names are exceeding 15 characters, 
then that is a sign that your code needs better organization.  
- Use of 1-character variables like `i,j,k` or `x` is acceptable only in loops
or in temporary variables in `lambda` functions.  
- Doc-strings should be concise and convey the most important information. 
In-line comments should be use sparingly when non-obvious logic is used. 

In general, when given the choice between retaining backward compatibility 
(coding for the past) or having the community update their code with the 
addition of a new library (coding for the future), I prefer the latter.

# Code of Conduct

Be courteous and respectful when providing (or receiving) feedback or 
criticism. Critique the code, not the person. 