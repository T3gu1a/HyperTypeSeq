# Hypergeometric-Type Sequences in Maple

**HyperTypeSeq** is a [Maple](https://www.maplesoft.com/) package to work with Hypergeometric-Type sequences. These are sequences whose general terms are linear combination of interlaced hypergeometric terms. The package provides:

- **HolonomiRE**: adapts $\texttt{HolonomicDE}$ from [FPS](https://www.mathematik.uni-kassel.de/~bteguia/FPS_webpage/FPS.htm) to search for a holonomic recurrence equation from an expression and a given bound. The syntax is
    $$\texttt{HolonomicRE}(\texttt{expr},\texttt{a}(\texttt{n}),\texttt{maxreorder}=\texttt{d},\texttt{reshift}=\texttt{t}),$$
where $\texttt{maxreorder}$ and $\texttt{reshift}$ are optional with default values $10$ and $1$, respectively. $\texttt{expr}$ is a term in $\texttt{n}$, and $\texttt{a}$ is the name of the unknown for the equation. $\texttt{maxreorder}$ is the maximum order of the holonomic recurrence equation sought, and $\texttt{reshift}$ is the minimal possible shift of $\texttt{a}(\texttt{n})$ in the recurrence equation sought.
    \item \texttt{REtoHTS}
- **REtoHTS**: 'decide' whether a given holonomic term is of hypergeometric type or not. The syntax is
    $$\texttt{REtoHTS}(\texttt{RE},\texttt{a}(\texttt{n}),\texttt{P}).$$
$\texttt{RE}$ is the holonomic recurrence equation and $\texttt{a}(\texttt{n})$ is the unknown term in it. $\texttt{P}$ is a procedure for computing values of the sequence at any index. $\texttt{P}$ can also be a list of initial values; however, the list must contain the values of the evaluations of $\texttt{expr}$ starting from $0$.
- **HTS**: writes a given expression into a hypergeometric-type normal form. The syntax is
    $$\texttt{HTS}(\texttt{expr},\texttt{n}),$$
with self-explanatory arguments from the previous commands. The argument $\texttt{maxreorder}$ is also optional for $\texttt{HTS}$.

## Installation

One can use **HyperTypeSeq** in Maple by putting the file HyperTypeSeq.mla in your working directory and include the lines
```
  > restart;

  > libname:=currentdir(), libname:

  > with(NLDE)
```
at the beginning of your Maple worksheet (session). To avoid putting these three lines in all worksheets, one can read the help page of the $\texttt{libname}$ command. It is not possible to generate this library from the source provided in this repository. Feel free to contact the author if you have any questions.

## Requirements and Dependencies

The package works best with Maple versions 2019 - 2021. There are some issues with the recent releases. The problem seems to come from the linear system solver $\texttt{SolveTools:-Linear}$.

## Author

- [Bertrand Teguia Tabuguia](https://bertrandteguia.com), University of Oxford
- licence: GNU General Public Licence v3.0.

## Examples

MapleWorksheet-HyperTypeSeq-examples.mw is a Maple session with some examples. The expected outputs are presented in MapleWorksheet-HyperTypeSeq-examples-outputs.mw

## References

1. [Hypergeometric-Type Sequences](https://). Bertrand Teguia Tabuguia. January 2024.


