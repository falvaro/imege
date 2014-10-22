IMEGE: Image-based Mathematical Expression Global Error
=======================================================
IMEGE is a measure for automatic evaluation of mathematical expression
recognition. This metric is computed from the LaTeX codification of
two mathematical expressions. Representation ambiguity is avoided by
rendering the image of given mathematical expressions and a normalized
global error value is computed. A detailed description can be found in

- Francisco Álvaro, Joan-Andreu Sánchez, José-Miguel Benedí. *"An
  image-based measure for evaluation of mathematical expression
  recognition"*. Iberian Conference on Pattern Recognition and Image
  Analysis (IbPRIA), 2013. Madeira, Portugal, June 5-7.
- Francisco Álvaro, Joan-Andreu Sánchez, and José-Miguel
  Benedí. *"IMEGE: Image-based Mathematical Expression Global
  Error"*. DSIC-PRHLT Technical Report, Universitat Politècnica de
  València, 2011.

License
-------
*IMEGE* is released under the [GNU General Public License version 3.0 (GPLv3)] [1]


Instructions
------------
This evaluation tool is made up of two parts, the BIDM algorithm and the IMEGE
script that computes the error.

The BIDM algorithm is written in C and uses the PGM and PBM utility libraries
by Jef Poskanzer. The code is already included, so for compiling it simply:

 1. Obtain the package using git:

        $ git clone https://github.com/falvaro/imege.git

    Or [download it as a zip file] [2]

 2. Go to the BIDM directory containing the source code.

 3. Compile the BIDM program

        $ make

As a result, you will have the executable file "*bidm*" ready to compute
the precision/recall between two images according to the description of the metric.
It has several options that are displayed when the program is executed without
arguments, but for calculating the IMEGE error it is not necessary to know how
BIDM works.

A bash script is provided as a wrapper in order to compute the
IMEGE error directly from LaTeX expressions. The images of the
expressions are generated to comply with format requirements of BIDM
(PGM format, depth 256) and BIDM parameters are properly set. This script
has the following requirements:

 - Packages *imagemagick* and *netpbm* for image treatment (commonly
   available in linux repositories).
 - Software [l2p] [3] to easily render the image of a given LaTeX expression.

Finally, the IMEGE error between the math expressions $x^2+1$ and $x2+1$
can be computed as

    $ ./IMEGE '$x^2+1$' '$x2+1$'

returning an 25.20% visual error.

To ensure mathematical representation uniqueness it is recommended to
add '*\displaystyle*' to every expression in order to avoid the layout
differences between $\sum_{i=1}^N x$ and $\displaystyle \sum_{i=1}^N x$.
Also, other special cases has to be taken into account, like $\sum_x$
versus $\Sigma_x$.


Citations
---------
If you use *IMEGE* for your research, please cite the following
reference:

<pre>
@INPROCEEDINGS{falvaro13b,
 author       = {Francisco \'Alvaro and Joan-Andreu Sánchez and José-Miguel Benedí},
 title        = {An image-based measure for evaluation of mathematical expression recognition},
 booktitle    = {6th Iberian Conference on Pattern Recognition and Image Analysis (IbPRIA), LNCS 7887},
 year         = {2013},
 pages        = {682-690},
 publisher    = {Springer},
}
</pre>



[1]: http://www.gnu.org/licenses/gpl-3.0.html
[2]: https://github.com/falvaro/imege/archive/master.zip
[3]: http://redsymbol.net/software/l2p
