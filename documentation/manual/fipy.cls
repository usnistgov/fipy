\NeedsTeXFormat{LaTeX2e} %[1995/12/01]
\ProvidesClass{fipy}
	     [2004/10/07 Document class (FiPy manual)]

% Change the options here to get a different set of basic options, but only
% if you have to.
%
\LoadClass[letterpaper]{book}

\def\tocPDFlevel{-1}
\DeclareOption{tocAsPDFpart}{\def\tocPDFlevel{-1}}
\DeclareOption{tocAsPDFchapter}{\def\tocPDFlevel{0}}
\ProcessOptions

% \usepackage[text={7in,9.333in}]{geometry}
\RequirePackage{geometry}

\RequirePackage{crop}

\RequirePackage{alltt, parskip, boxedminipage} % fancyheadings,  
\RequirePackage{multirow, longtable, makeidx, amssymb} 

\RequirePackage[other]{tocbibind}
\tocotherhead{part}
\setlength{\columnsep}{35pt}

\RequirePackage{amsmath} 

\RequirePackage{minitoc}
% \setcounter{tocdepth}{0}
\setcounter{parttocdepth}{1}

\RequirePackage[ascii]{inputenc}

% \usepackage{fullpage}
%\usepackage[headings]{fullpage}
\makeindex
\RequirePackage[usenames]{color}
\definecolor{darkblue}{rgb}{0,0.05,0.35}
\definecolor{redish}{rgb}{0.894,0.122,0.122}
\definecolor{bluish}{rgb}{0.216,0.188,0.533}
% \usepackage[dvips, pagebackref, pdftitle={FiPy}, pdfcreator={epydoc 2.1}, bookmarks=true, bookmarksopen=false, pdfpagemode=UseOutlines, colorlinks=true, linkcolor=black, anchorcolor=black, citecolor=black, filecolor=black, menucolor=black, pagecolor=black, urlcolor=darkblue]{hyperref}
\RequirePackage[urlcolor=blue,linkcolor=blue,bookmarksopen,bookmarksopenlevel=2,pdftex,pagebackref,pdftitle={FiPy}, pdfcreator={epydoc 2.1}, bookmarks=true, bookmarksopen=false, pdfpagemode=UseOutlines,colorlinks=true,plainpages=false,pdfpagelabels]{hyperref}
\RequirePackage{graphicx}
% \usepackage{memhfixc}
% \usepackage{nameref}

\graphicspath{{../figures/},{images/}}

% for reStructuredText
\RequirePackage{shortvrb}
\RequirePackage{tabularx}
\setlength{\extrarowheight}{2pt}
\newlength{\admonitionwidth}
\setlength{\admonitionwidth}{0.9\textwidth}
\newlength{\docinfowidth}
\setlength{\docinfowidth}{0.9\textwidth}
\newlength{\locallinewidth}
\newcommand{\optionlistlabel}[1]{\bf #1 \hfill}
\newenvironment{optionlist}[1]
{\begin{list}{}
  {\setlength{\labelwidth}{#1}
   \setlength{\rightmargin}{1cm}
   \setlength{\leftmargin}{\rightmargin}
   \addtolength{\leftmargin}{\labelwidth}
   \addtolength{\leftmargin}{\labelsep}
   \renewcommand{\makelabel}{\optionlistlabel}}
}{\end{list}}
% begin: floats for footnotes tweaking.
\setlength{\floatsep}{0.5em}
\setlength{\textfloatsep}{\fill}
\addtolength{\textfloatsep}{3em}
\renewcommand{\textfraction}{0.5}
\renewcommand{\topfraction}{0.5}
\renewcommand{\bottomfraction}{0.5}
\setcounter{totalnumber}{50}
\setcounter{topnumber}{50}
\setcounter{bottomnumber}{50}
% end floats for footnotes
% some commands, that could be overwritten in the style file.
\newcommand{\rubric}[1]{\subsection*{~\hfill {\it #1} \hfill ~}}
\newcommand{\titlereference}[1]{\textsl{#1}}


\newcommand{\logo}{\rotatebox{4}{\textcolor{redish}{\( \varphi \)}}\kern-.70em\raisebox{-.15em}{\textcolor{bluish}{\( \pi\)}}}
% \newcommand{\logoToo}{\raisebox{-.15em}{\textcolor{bluish}{\(\pi\)}}\kern-.64em\rotatebox{4}{\textcolor{redish}{\( \varphi \)}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% \includeonly{installation}


\renewcommand*\l@section{\@dottedtocline{1}{1.5em}{3.3em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{3.8em}{4.2em}}
\renewcommand*\l@subsubsection{\@dottedtocline{3}{7.0em}{5.1em}}
\renewcommand*\l@paragraph{\@dottedtocline{4}{10em}{6em}}
\renewcommand*\l@subparagraph{\@dottedtocline{5}{12em}{7em}}

% \usepackage{layout}

\RequirePackage{tocvsec2}

\input{version}

\setlength{\fboxrule}{2\fboxrule}
\newlength{\BCL} % base class length, for base trees.

% \pagestyle{Ruled}
\renewcommand{\sectionmark}[1]{\markboth{#1}{}}
\renewcommand{\subsectionmark}[1]{\markright{#1}}

\newenvironment{Ventry}[1]%
  {\begin{list}{}{%
    \renewcommand{\makelabel}[1]{\texttt{##1:}\hfil}%
    \settowidth{\labelwidth}{\texttt{#1:}}%
    \setlength{\leftmargin}{\labelsep}%
    \addtolength{\leftmargin}{\labelwidth}}}%
  {\end{list}}

\newsavebox{\NISTbox}
\sbox{\NISTbox}{\includegraphics[trim=5 2 5 5,scale=1.5]{NIST_right2line}}
%\newlength{\photoheight}
%\settoheight{\photoheight}{\usebox{\photobox}}
%\newlength{\photowidth}
%\settowidth{\photowidth}{\usebox{\photobox}}

\newenvironment{changemargin}[2]{%
 \begin{list}{}{%
  \setlength{\topsep}{0pt}%
  \setlength{\leftmargin}{#1}%
  \setlength{\rightmargin}{#2}%
  \setlength{\listparindent}{\parindent}%
  \setlength{\itemindent}{\parindent}%
  \setlength{\parsep}{\parskip}%
 }%
\item[]}{\end{list}}

\renewcommand\tableofcontents{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\contentsname\pdfbookmark[\tocPDFlevel]{Contents}{contents}
	\@mkboth{%
	   \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
    }

\newcommand*{\ps@mytitlepage}{%
  \let\@oddhead\@empty
  \renewcommand*{\@oddfoot}{%
  \rule{6.4in}{0in}\raisebox{-0.6in}{\parbox[b][0pt]{0pt}{\makebox[0pt][r]{\usebox{\NISTbox}}}}}}%\vspace*{-0.5in}}

\def\@makeauthor{\begingroup
\renewcommand{\and}{\\}
\begin{tabular}[t]{r@{}}%
    \@author
\end{tabular}
\endgroup
}

\def\@makeaffiliation{\begingroup
\renewcommand{\and}{\\}
\textit{
\begin{tabular}[t]{r@{}}%
    \@affiliation
\end{tabular}}
\endgroup
}

\def\subtitle#1{\gdef\@subtitle{#1}}
\def\@subtitle{}
\def\affiliation#1{\begingroup
\renewcommand{\and}{\\}
\gdef\@affiliation{#1}
\endgroup
}
\def\@affiliation{}

\renewcommand{\maketitle}{%
\begin{titlepage}%
\thispagestyle{mytitlepage}
\begin{changemargin}{-1in}{-0.5in}
\let\footnotesize\small
\let\footnoterule\relax
\let \footnote \thanks
\begin{flushright}%
  \vspace*{-0.6in}
  \scalebox{10}{\logo}\par
  \vskip 3em%
  {\Huge \@title \par}%
  {\huge \@subtitle \par}%
  \vskip 3em%
  {\large
   \lineskip .75em%
   \@makeauthor \par
   \@makeaffiliation \par }%
    \vskip 1.5em%
  {\large \@date \par}%       % Set date in \large size.
      \vskip 1.5em%
  {\large Version~\Version \par}
  \vspace*{\fill}%
\end{flushright}%
\end{changemargin}%
\end{titlepage}%
}

\newcommand{\@Code}[1]{\texttt{#1}}
\newcommand{\Class}[1]{\@Code{#1}\index{#1@\@Code{#1} class}}

% commands for various weblinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\@Software}[1]{\textsf{#1}}
\newcommand{\@Org}[1]{\textsf{#1}}

\newcommand{\FiPy}{\@Software{FiPy}}

\newcommand{\@OrgURL}[2]{\href{#2}{\@Org{#1}}\index{#1@\@Org{#1}}}

\newcommand{\@SoftwareURL}[2]{\href{#2}{\@Software{#1}}\index{#1@\@Software{#1}}}

\newcommand{\SciPy}{\@SoftwareURL{SciPy}{http://www.scipy.org/}}

\newcommand{\PyRex}{\@SoftwareURL{PyRex}{http://www.cosc.canterbury.ac.nz/~greg/python/Pyrex/}}

\newcommand{\Python}{\@SoftwareURL{Python}{http://www.python.org/}}

\newcommand{\Numeric}{\@SoftwareURL{Numeric}{http://www.numpy.org/}}

\newcommand{\PyX}{\@SoftwareURL{PyX}{http://pyx.sourceforge.net/}}

\newcommand{\PyGist}{\@SoftwareURL{PyGist}{http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/pygist.html}}

\newcommand{\PySparse}{\@SoftwareURL{PySparse}{http://people.web.psi.ch/geus/pyfemax/pysparse.html}}

\newcommand{\ScientificPython}{\@SoftwareURL{Scientific Python}{http://starship.python.net/~hinsen/ScientificPython/}}

\newcommand{\MacOSX}{\@SoftwareURL{Mac OS X}{http://www.apple.com/macosx/}}

\newcommand{\Gmsh}{\@SoftwareURL{Gmsh}{http://www.geuz.org/gmsh/}}

\newcommand{\weave}{\@SoftwareURL{weave}{http://www.scipy.org/documentation/weave/}}

\newcommand{\blitz}{\@SoftwareURL{blitz}{http://www.scipy.org/documentation/weave/}}

\newcommand{\doctest}{\@SoftwareURL{doctest}{http://www.python.org/doc/current/lib/module-doctest.html}}

\newcommand{\unittest}{\@SoftwareURL{unittest}{http://www.python.org/doc/current/lib/module-unittest.html}}

\newcommand{\NIST}{\@OrgURL{NIST}{http://www.nist.gov/}}

\newcommand{\CTCMS}{\@OrgURL{http://www.ctcms.nist.gov/}{\textsf{CTCMS}}}

\newcommand{\MSEL}{\@OrgURL{http://www.msel.nist.gov/}{\textsf{MSEL}}}


