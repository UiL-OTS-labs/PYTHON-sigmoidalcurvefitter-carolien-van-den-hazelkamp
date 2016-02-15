#sigcurvefitter
Provides a least-square curve fitting implementation coded in python for a sigmoidal curve that fits the sigmoidal curve equation:

`d(time) = lambda*(1-math.exp(-beta*(time-delta)))`
Where if time <= delta 
* En dan is lambda de asymptote, delta de intercept met de x-as (waar accuracy van kans afwijkt) en beta de steilheid van de curve. *
If time > delta, d(time) is 0. 


# Dependencies
The script relies on python and several libraries. Likely libraries and the depencies you'll need to install are matplotlib, numpy and scipy. The best suggestion is to download a total package. I suggest Canopy (Express free)
 `https://store.enthought.com/`.

# Usage
 `python sigcurvefitter.py arg1`
Where arg1 is an optional CSV filename (See *exampledata.csv* for the format of the input document)

# Authors
Internals are inspired by 
 `http://people.duke.edu/~ccc14/pcfb/analysis.html`

Coded by **Jan de Mooij** and fine-tuned by **Chris van Run**
