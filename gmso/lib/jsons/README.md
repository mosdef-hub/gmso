# GMSO Potential Templates JSON Files

This directory contains json representations of GMSO potential 
templates. 

A potential template object requires the following attributes:

1. `name`: Name of the potential template
2. `expression`: expression for the potential template
3. `independent_variables`: Independent variables for the template 


## For Developers
To add/contribute a new `PotentialTemplate` there are two ways:

1. **Add JSON File**: The example below shows the content of an example json file(OPLSTorsionPotential):
```json
{
  "name": "OPLSTorsionPotential",
  "expression": "0.5 * k0 + 0.5 * k1 * (1 + cos(phi)) + 0.5 * k2 * (1 - cos(2*phi)) + 0.5 * k3 * (1 + cos(3*phi)) + 0.5 * k4 * (1 - cos(4*phi))",
  "independent_variables": "phi"
}
```
If you want to add a template named `ABCTemplate`. Create a file called
`ABCTemplate.json` in this directory. Using the json reference above, 
name the template `ABCTemplate` and change the `expression` and `independent_variables`.

If there are more than one  `independent_variables` in your template, simply separate them by a comma(,).

2. **Use a python script**:
Let's say you want to create a new PotentialTemplate. There's a static method 
in `PotentialTemplates` that let's you add a new one.

Use the following script to generate potential template(in your dev branch):

**Warning: Do not use this script on conda/pip installed version or master branch**

```python
from gmso.lib.potential_templates import PotentialTemplates
templates = PotentialTemplates()
abcTemplate_dict = {
            'name': 'ABCTemplate', 
            'expression': 'expression', 
            'independent_variables': 'variables'
        }       
templates.save_potential_template('ABCTemplate', abcTemplate_dict, True)
```