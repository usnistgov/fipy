\usepackage[text={7in,9.333in}]{geometry}

\usepackage{crop}

\usepackage{alltt, parskip, boxedminipage} % fancyheadings,  
\usepackage{multirow, longtable, makeidx, tocbibind, amssymb} 
\usepackage{amsmath} 

\usepackage{minitoc}
% \setcounter{tocdepth}{0}
\setcounter{parttocdepth}{1}

\usepackage[ascii]{inputenc}

%\usepackage{multirow, longtable, tocbibind, amssymb} % makeidx, 
% \usepackage{fullpage}
%\usepackage[headings]{fullpage}
\makeindex
\usepackage[usenames]{color}
\definecolor{darkblue}{rgb}{0,0.05,0.35}
\definecolor{redish}{rgb}{0.894,0.122,0.122}
\definecolor{bluish}{rgb}{0.216,0.188,0.533}
% \usepackage[dvips, pagebackref, pdftitle={FiPy}, pdfcreator={epydoc 2.1}, bookmarks=true, bookmarksopen=false, pdfpagemode=UseOutlines, colorlinks=true, linkcolor=black, anchorcolor=black, citecolor=black, filecolor=black, menucolor=black, pagecolor=black, urlcolor=darkblue]{hyperref}
\usepackage[urlcolor=blue,linkcolor=blue,bookmarksopen,bookmarksopenlevel=2,pdftex,pagebackref,pdftitle={FiPy}, pdfcreator={epydoc 2.1}, bookmarks=true, bookmarksopen=false, pdfpagemode=UseOutlines,colorlinks=true,plainpages=false,pdfpagelabels]{hyperref}
\usepackage{graphicx}
% \usepackage{memhfixc}
% \usepackage{nameref}

\graphicspath{{../figures/}}

% for reStructuredText
\usepackage{shortvrb}
\usepackage{tabularx}
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

% commands for various weblinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\FiPy}{\textsf{FiPy}}

\newcommand{\SciPy}{\href{http://www.scipy.org/}{\textsf{SciPy}}}

\newcommand{\PyRex}{\href{http://www.cosc.canterbury.ac.nz/~greg/python/Pyrex/}{\textsf{PyRex}}}

\newcommand{\Python}{\href{http://www.python.org/}{\textsf{Python}}}

\newcommand{\Numeric}{\href{http://www.numpy.org/}{\textsf{Numeric}}}

\newcommand{\NIST}{\href{http://www.nist.gov/}{\textsf{NIST}}}

\newcommand{\CTCMS}{\href{http://www.ctcms.nist.gov/}{\textsf{CTCMS}}}

\newcommand{\MSEL}{\href{http://www.msel.nist.gov/}{\textsf{MSEL}}}

\newcommand{\PyX}{\href{http://pyx.sourceforge.net/}{\textsf{PyX}}}

\newcommand{\PyGist}{\href{http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/python/pygist.html}{\textsf{PyGist}}}

\newcommand{\MacOSX}{\href{http://www.apple.com/macosx/}{\textsf{Mac OS X}}}

\newcommand{\Gmsh}{\href{http://www.geuz.org/gmsh/}{\textsf{Gmsh}}}

\newcommand{\weave}{\href{http://www.scipy.org/documentation/weave/}{\textsf{weave}}}

\newcommand{\blitz}{\href{http://www.scipy.org/documentation/weave/}{\textsf{blitz}}}

\newcommand{\doctest}{\href{http://www.python.org/doc/current/lib/module-doctest.html}{\textsf{doctest}}}

\newcommand{\unittest}{\href{http://www.python.org/doc/current/lib/module-unittest.html}{\textsf{unittest}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% \includeonly{installation}

\makeatletter
\renewcommand*\l@section{\@dottedtocline{1}{1.5em}{3.3em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{3.8em}{4.2em}}
\renewcommand*\l@subsubsection{\@dottedtocline{3}{7.0em}{5.1em}}
\renewcommand*\l@paragraph{\@dottedtocline{4}{10em}{6em}}
\renewcommand*\l@subparagraph{\@dottedtocline{5}{12em}{7em}}

\renewcommand\tableofcontents{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\contentsname\pdfbookmark[-1]{Contents}{contents}
	\@mkboth{%
	   \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
    }

\renewcommand\maketitle{\begin{titlepage}%
\let\footnotesize\small
\let\footnoterule\relax
\let \footnote \thanks
\null\vfil
\vskip 20\p@
\begin{flushright}%
  {\Huge \@title \par}%
  \vskip 3em%
  {\large
   \lineskip .75em%
    \begin{tabular}[t]{r@{}}%
      \@author
    \end{tabular}\par}%
    \vskip 1.5em%
  {\large \@date \par}%       % Set date in \large size.
\end{flushright}\par
\@thanks
\vfil\null
\end{titlepage}%
\setcounter{footnote}{0}%
\global\let\thanks\relax
\global\let\maketitle\relax
\global\let\@thanks\@empty
\global\let\@author\@empty
\global\let\@date\@empty
\global\let\@title\@empty
\global\let\title\relax
\global\let\author\relax
\global\let\date\relax
\global\let\and\relax
}

\makeatother

% \usepackage{layout}


\usepackage{tocvsec2}

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

\makeatletter
\newcommand*{\ps@mytitlepage}{%
  \let\@oddhead\@empty
  \renewcommand*{\@oddfoot}{%
  \hfill\parbox[b][0pt]{0pt}{\makebox[0pt][r]{\usebox{\NISTbox}}}}}%\vspace*{-0.5in}}

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
\let\footnotesize\small
\let\footnoterule\relax
\let \footnote \thanks
\begin{flushright}%
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
\end{titlepage}%
}

\makeatother
