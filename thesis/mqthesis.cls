% mqthesis.cls - Macquarie University PhD and MRes Thesis LaTeX Class File.
%
% Copyright (C) 2004 Paul Cochrane
% Copyright (C) 2008-2015 Alexei Gilchrist
% 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%   
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%    
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mqthesis}[2013/10/14 Macquarie University physics thesis LaTeX class]

% MQ2015 Closest match for Georgia with math support:
%\usepackage{newtxtext}
% newtxmath is loaded below, after amsmath


% ifthenelse for if loops
\RequirePackage{ifthen}

% is the thesis for Phd, MRes or Hons?  
% define the default first
\newcommand{\degreetext}{for the degree of Doctor of Philosophy}
\DeclareOption{phd}{\renewcommand{\degreetext}%
{for the degree of Doctor of Philosophy}}
\DeclareOption{mres}{\renewcommand{\degreetext}%
{for the degree of Master of Science (Mathematics and Statistics)}}
\DeclareOption{hons}{\renewcommand{\degreetext}%
{in partial fulfilment of the degree of Bachelor of Science with Honours}}

% print a copyright page on the back of the title page?
\newboolean{copyrightpage}
\setboolean{copyrightpage}{false}
\DeclareOption{copyrightpage}{\setboolean{copyrightpage}{true}}

% print ``Examiner's copy'' on the title page?
\newboolean{examinerscopy}
\setboolean{examinerscopy}{false}
\DeclareOption{examinerscopy}{\setboolean{examinerscopy}{true}}


% print grey chapter numbers instead of MQ corporate color "sand"?
\newboolean{greychapternumbers}
\setboolean{greychapternumbers}{false}
\DeclareOption{greychapternumbers}{\setboolean{greychapternumbers}{true}}

% finished processing options
\ProcessOptions\relax

% load the book class using a4paper
\LoadClass[12pt, a4paper]{book}

% MQ specification 
% "Margins should be not less than 3.5cm on the binding edge, 
% 1.5cm on the opposite edge and 1.5cm at the top and bottom 
% to allow for binding and trimming.""
\usepackage[top=3.0cm, bottom=3.0cm, inner=3.0cm, outer=3.0cm, includehead]{geometry}


\brokenpenalty=10000   % dunno what this does, maybe handy

% this stops one figure taking up a whole page and lets more text onto
% the one page when a figure exists
\renewcommand\floatpagefraction{0.8} %   Default = 0.5

% macros for defining the department, various dates and the title and
% author of the thesis.
\def\department#1{\gdef\@department{#1}}
\def\submitdate#1{\gdef\@submitdate{#1}}
\def\copyrightyear#1{\gdef\@copyrightyear{#1}} % \author, \title in report
\def\@title{Logic-Based Benders Decomposition applied to the Setup Assembly Line Balancing and Scheduling Problem}
\def\@author{Kenneth D. Young}
\def\@supervisors{Dr. Alysson M. Costa}
\def\@submitdate{\ifcase\the\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space \number\the\year}
\def\@copyrightyear{\number\the\year}

%%% load the required packages
% fancyhdr for nice, fancy headings
\RequirePackage{fancyhdr}
% xspace so that spaces after commands are handled correctly
\RequirePackage{xspace}

% This constructs the MQ title page
\def\titlep{%
    \thispagestyle{empty}%
    \null\vskip.5in%
    \begin{center}
            \begin{doublespacing}
                {\huge\expandafter{\textbf{\@title}}} \\[2.5cm]
            \end{doublespacing}
            By \\[1.0cm]
            Kenneth D. Young\\[0.75cm]
            Supervisor: \text{\@supervisors}\\
    \end{center}
    \vfill
    \begin{center}
        A thesis submitted to the University of Melbourne\\
        \degreetext{}\\
        \expandafter{\@department}\\
        \@submitdate\\[3cm]
        %\includegraphics[width=0.3 \columnwidth]{MQlogo}\\
        %\includegraphics[width=0.3 \columnwidth]{UnimelbLogo}\\
    \end{center}
    \ifthenelse{\boolean{examinerscopy}}{%
        \begin{center}
            %\vskip1in
            {\LARGE {Examiner's Copy}}
        \end{center}%
        }{}
    \vfill
    \newpage}

% this defines the signature for the signature page where the standard
% declaration is put
\long\def\signature#1{%
\begin{center}
\begin{minipage}{4in}
\parindent=0pt
Except where acknowledged in the customary manner, the material
presented in this thesis is, to the best of my knowledge, original and
has not been submitted in whole or part for a degree in any university.
\par
\vspace{.5in}
\hbox to 4in{\hfil\shortstack{\vrule width 3in height 0.4pt\\#1}}
\end{minipage}
\end{center}}

% this is the actual signature page
\def\signaturepage{%
	\vskip.5in
	\signature{\@author}}

% now put in the titlepage
% (bit of a misnomer really, as it is really the first three or four pages.
\def\titlepage{%
    \pagenumbering{roman}  % make sure that front matter is numbered Roman
    \titlep
    \ifthenelse{\boolean{copyrightpage}}{%
        \newpage
        \vspace*{\fill}
        \begin{center}
            \copyright{} \@author, \@copyrightyear.
        \end{center}
        \begin{center}
            Typeset in \LaTeXe{}.
        \end{center}
        \normalsize
        \newpage
        \signaturepage}{}
}


\usepackage[pdftex,usenames,dvipsnames]{color}

\usepackage[grey,times]{quotchap}

% set the pagestyle to look good
\pagestyle{fancy}

%%%%% Fancyhdr stuff
% give the header a bit more room, otherwise LaTeX will spew on each page
\addtolength{\headheight}{2.5pt}
% define how headers are marked, for details, see fancyhdr docs
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}

% define where sections, chapters and pagenumbers are put
% see fancyhdr docs for details
% the \nouppercase stops book.cls making the contents, bibliography
% and index headers from being all in uppercase.
% The options used here are essentially that in Lamport's book, but
% with small caps for the headings.
\fancyhf{}
\fancyhead[LE,RO]{\nouppercase{\thepage}}
\fancyhead[LO]{\nouppercase{\rightmark}}
\fancyhead[RE]{\nouppercase{\leftmark}}

% \chapter issues a \thispagestyle{plain} .. override plain
% Redefine the plain page style
\fancypagestyle{plain}{%
  \fancyhf{}%
   \renewcommand{\headrulewidth}{0pt}% Line at the header invisible
}

%%% other settings required for a thesis
% It's a references section, not a bibliography, hence redefine
% \bibname i.e. change ``Bibliography'' to ``References''
\renewcommand*{\bibname}{References}

% use 1.5 line spacing for the examiners copy
% single line spacing for the final copy
\ifthenelse{\boolean{examinerscopy}}{%
\renewcommand{\baselinestretch}{1.5}}%
{\renewcommand{\baselinestretch}{1}}

% spell things correctly
\newenvironment{centre}{\begin{center}}{\end{center}}
\newenvironment{itemise}{\begin{itemize}}{\end{itemize}}

\usepackage{makeidx}

% MQ2015  Sand=R214 G210 B196  
\ifthenelse{\boolean{greychapternumbers}}{}{%
\definecolor{chaptergrey}{rgb}{0.839,0.824,0.769}%
}

%%%%% set up the bibliography style
\bibliographystyle{myplainnat}  % mqthesis bibliography style file, made
			                  % with makebst



%----------------------------------------------%
% Packages %
%----------------------------------------------%

%%%%% optional packages
\usepackage[square,comma,numbers,sort&compress]{natbib}
		% this is the natural sciences bibliography citation
		% style package.  The options here give citations in
		% the text as numbers in square brackets, separated by
		% commas, citations sorted and consecutive citations
		% compressed 
		% output example: [1,4,12-15]

\usepackage[]{tocbibind}  
				% allows the table of contents, bibliography
				% and index to be added to the table of
				% contents if desired, the option used
				% here specifies that the table of
				% contents is not to be added.
				% tocbibind needs to be after natbib
				% otherwise bits of it get trampled.

\usepackage{color}
\usepackage{graphicx} % standard graphics package for inclusion of
              % images and eps files into LaTeX document

\usepackage{amsmath,amsfonts,amssymb,amsthm} % this is handy for mathematicians and physicists
			      % see http://www.ams.org/tex/amslatex.html
\usepackage{latexsym,epsfig,url,fancyhdr,verbatim}
\usepackage{float}
\usepackage{epstopdf}
\usepackage{IEEEtrantools}
\usepackage{multicol,multirow}
\usepackage{authblk}
\usepackage{chngpage}
\usepackage{mathrsfs}
\usepackage{listings}
\usepackage{booktabs}
\usepackage{setspace}
\usepackage{xifthen}
\usepackage{algpseudocode}
\usepackage{algorithm}
\usepackage{algorithmicx}
\usepackage{mdframed}
%\usepackage{xr}


%----------------------------------------------%
% Misc. %
%----------------------------------------------%
\theoremstyle{definition}
\newmdtheoremenv[
  hidealllines=true,
  leftline=true,
  innerleftmargin=10pt,
  innerrightmargin=10pt,
  innertopmargin=0pt,
]{example}{Example}%{\itshape}{\rmfamily}
\newmdtheoremenv[
  hidealllines=true,
  innerleftmargin=10pt,
  innerrightmargin=10pt,
  innertopmargin=0pt,
]{theorem}{Theorem}
\newmdtheoremenv[
  hidealllines=true,
  innerleftmargin=10pt,
  innerrightmargin=10pt,
  innertopmargin=0pt,
]{lemma}{Lemma}


%----------------------------------------------%
% Author's comments %
%----------------------------------------------%
\setlength\marginparwidth{2cm}
\newcommand{\ky}[1]%
{\marginpar{\textcolor{blue}{\textsc{Kenneth}}}\textcolor{blue}{[[#1]]}}
\newcommand{\todo}[1]%
{\marginpar{\textcolor{red}{\textsc{ToDo}}}\textcolor{red}{[[#1]]}}

%----------------------------------------------%
% Citation Commands %
%----------------------------------------------%
\newcommand{\authcite}[1]{\citeauthor{#1} \citeyear{#1}}
\newcommand{\authciteb}[1]{\citeauthor{#1} (\citeyear{#1})}

%----------------------------------------------%
% Algorithmic Package %
%----------------------------------------------%
\newcommand*\Let[2]{\State #1 $\gets$ #2}
\algrenewcommand\algorithmicrequire{\textbf{Precondition:}}
\algrenewcommand\algorithmicensure{\textbf{Postcondition:}}

\algtext*{EndWhile}% Remove "end while" text
\algtext*{EndIf}% Remove "end if" text
\algtext*{EndFor}% Remove "end for" text

%----------------------------------------------%
% Notation Commands %
%----------------------------------------------%
\newcommand{\aldif}{\texttt{alldifferent}\xspace}
\newcommand{\cumu}{\texttt{cumulative}\xspace}
\newcommand{\disj}{\texttt{disjunctive}\xspace}

\newcommand{\eg}{\emph{e.g.},\xspace}
\newcommand{\ie}{\emph{i.e.},\xspace}
\newcommand{\cf}{\emph{cf.}\xspace}
\newcommand{\ap}{\emph{a priori},\xspace}

\newcommand{\chuffed}{\texttt{Chuffed}\xspace}
\newcommand{\gurobi}{Gurobi\xspace}

\newcommand{\CP}{CP\xspace}
\newcommand{\LCG}{LCG\xspace}
\newcommand{\MSPSP}{MSPSP\xspace}
\newcommand{\RCPSP}{RCPSP\xspace}
\newcommand{\SAT}{SAT\xspace}

\newcommand{\intser}{\texttt{int\_search}\xspace}
\newcommand{\bolser}{\texttt{bool\_search}\xspace}
\newcommand{\seqser}{\texttt{seq\_search}\xspace}
\newcommand{\priser}{\texttt{priority\_search}\xspace}

\newcommand{\bab}{B\&B\xspace}

% Problems
\newcommand{\albp}[1]{\ifthenelse{\isempty{#1}}{ALBP\xspace}{ALBP-#1\xspace}}
\newcommand{\sab}[1]{\ifthenelse{\isempty{#1}}{SALBP\xspace}{SALBP-#1\xspace}}
\newcommand{\sua}[1]{\ifthenelse{\isempty{#1}}{SUALBSP\xspace}{SUALBSP-#1\xspace}}
\newcommand{\gab}[1]{\ifthenelse{\isempty{#1}}{GALBSP\xspace}{GALBSP-#1\xspace}}
\newcommand{\suwa}[1]{\ifthenelse{\isempty{#1}}{SUALWABSP\xspace}{SUALWABSP-#1\xspace}}
\newcommand{\suwi}[1]{\ifthenelse{\isempty{#1}}{SUALWIBSP\xspace}{SUALWIBSP-#1\xspace}}

% MIPs
\newcommand{\fsbf}[1]{\ifthenelse{\isempty{#1}}{FSBF\xspace}{FSBF-#1\xspace}}
\newcommand{\ssbf}[1]{\ifthenelse{\isempty{#1}}{SSBF\xspace}{SSBF-#1\xspace}}
\newcommand{\scbf}[1]{\ifthenelse{\isempty{#1}}{SCBF\xspace}{SCBF-#1\xspace}}

% Benders problems
\newcommand{\rmp}[1]{\ifthenelse{\isempty{#1}}{RMP\xspace}{RMP($#1$)\xspace}}
\newcommand{\spmip}[1]{\ifthenelse{\isempty{#1}}{SP-MIP\xspace}{SP-MIP($#1$)\xspace}}
\newcommand{\spcp}[1]{\ifthenelse{\isempty{#1}}{SP-CP\xspace}{SP-CP($#1$)\xspace}}
\newcommand{\sptsp}[1]{\ifthenelse{\isempty{#1}}{SP-TSP\xspace}{SP-TSP($#1$)\xspace}}

\newcommand{\grasp}{GRASP\xspace}


\newcommand{\ol}[1]{\overline{#1}\xspace}
\newcommand{\ul}[1]{\underline{#1}\xspace}

%----------------------------------------------%
% IEEEeqnarray settings %
%----------------------------------------------%
\interdisplaylinepenalty=1000


%----------------------------------------------%
% Define Lengths, Spacings, etc. %
%----------------------------------------------%
\setstretch{1.0}
\setlength{\parskip}{4pt}

% Define spacing between IEEEeqnarray lines
\newcommand{\eqnv}{4pt}
\newcommand{\smolEqnv}{0pt}


%----------------------------------------------%
% Subsubsection counters
%----------------------------------------------%
\setcounter{tocdepth}{4}
\setcounter{secnumdepth}{4}

%----------------------------------------------%
% Syntax highlighting for MiniZinc in listings %
%----------------------------------------------%

\definecolor{ForestGreen}{RGB}{34,106,46}
\definecolor{lightgray}{rgb}{0.97, 0.97, 0.97}

\lstdefinelanguage{minizinc}{
    morekeywords={
        %% MiniZinc keywords
        %%
        ann, annotation, any, array, assert,
        bool,
        constraint,
        else, elseif, endif, enum, exists,
        float, forall, function,
        if, in, include, int,
        list,
        minimize, maximize,
        of, op, output,
        par, predicate,
        record,
        set, solve, string,
        test, then, tuple, type,
        var,
        where,
        %% MiniZinc functions
        %%
        abort, abs, acosh, array_intersect, array_union,
        array1d, array2d, array3d, array4d, array5d, array6d, asin, assert, atan,
        bool2int,
        card, ceil, combinator, concat, cos, cosh,
        dom, dom_array, dom_size, dominance,
        exp,
        fix, floor,
        index_set, index_set_1of2, index_set_2of2, index_set_1of3, index_set_2of3, index_set_3of3,
        int2float, is_fixed,
        join,
        lb, lb_array, length, let, ln, log, log2, log10,
        min, max,
        pow, product,
        round,
        set2array, show, show_int, show_float, sin, sinh, sqrt, sum,
        tan, tanh, trace,
        ub, and ub_array,
        %% Search keywords
        %%
        bool_search, int_search, seq_search, priority_search,
        %% MiniSearch keywords
        %%
        minisearch, search, while, repeat, next, commit, print, post, sol, scope, time_limit, break, fail
    },
    sensitive=true, % are the keywords case sensitive
    morecomment=[l][\em\color{ForestGreen}]{\%},
    %morecomment=[s]{/*}{*/},
    morestring=[b]",
}

%% Settings for listings
%%
\lstset{ %
    backgroundcolor=\color{lightgray},  % choose the background color; you must add
                                        % \usepackage{color} or \usepackage{xcolor}
    basicstyle=\scriptsize\ttfamily,    % the size of the fonts that are used for the code
    belowskip=-2em,
    breakatwhitespace=false,            % sets if automatic breaks should only happen at whitespace
    breaklines=true,                    % sets automatic line breaking
    captionpos=b,                       % sets the caption-position to bottom
    commentstyle=\color{ForestGreen},   % comment style
    %deletekeywords={...},              % if you want to delete keywords from the given language
    escapeinside={\%*}{*)},             % if you want to add LaTeX within your code
    extendedchars=true,                 % lets you use non-ASCII characters; for 8-bits
                                        % encodings only, does not work with UTF-8
    frame=single,                       % adds a frame around the code
    keepspaces=true,                    % keeps spaces in text, useful for keeping indentation
                                        % of code (possibly needs columns=flexible)
    keywordstyle=\bfseries\color{blue}, % keyword style
    language=minizinc,                  % the language of the code
    %morekeywords={*,...},              % if you want to add more keywords to the set
    numbers=none,                       % where to put the line-numbers; possible values are (none, left, right)
    %numbersep=5pt,                     % how far the line-numbers are from the code
    %numberstyle=\tiny\color{Gray},     % the style that is used for the line-numbers
    rulecolor=\color{black},            % if not set, the frame-color may be changed
                                        % on line-breaks within not-black text (e.g. comments (green here))
    showspaces=false,                   % show spaces everywhere adding particular
                                        % underscores; it overrides 'showstringspaces'
    showstringspaces=false,             % underline spaces within strings only
    showtabs=false,                     % show tabs within strings adding particular underscores
    %stepnumber=1,                      % the step between two line-numbers. If it's 1, each line will be numbered
    stringstyle=\color{Red},            % string literal style
    tabsize=2,                          % sets default tabsize to 2 spaces
    title=\lstname                      % show the filename of files included with \lstinputlisting;
                                        % also try caption instead of title
}

%\def\mzninline{\lstinline[basicstyle=\ttfamily,annotationstyle=\normalfont]}
\def\mzninline{\verb}



% load this after amsmath
%\usepackage{newtxmath}

%\usepackage[pdfusetitle]{hyperref}
%\hypersetup{
%    colorlinks=false,  % false: boxed links; true: colored links
%    linkcolor=black,   % color of internal links 
%    citecolor=black,   % color of links to bibliography
%    filecolor=black,   % color of file links
%    urlcolor=black,    % color of external links
%    }
%\urlstyle{same}


