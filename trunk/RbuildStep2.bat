PATH = C:\programs\R\R-3.2.0\bin;C:\programs\R\R-3.2.0\bin\x64;C:\programs\R\R-3.2.0\share\texmf;C:\programs\R\library;C:\Rtools\bin;C:\Rtools\gcc-4.6.3\bin;C:\Rtools\gcc-4.6.3\bin64;C:\Rtools\gcc-4.6.3\include;C:\programs\MiKTeX 2.9\miktex\bin;C:\programs\Ghostgum64\gsview;C:\programs\gs64\gs9.06\bin;C:\Program Files (x86)\Adobe\Acrobat 10.0\Acrobat;C:\Program Files\HTML Help Workshop;C:\WINDOWS;%HOMEDRIVE%%HOMEPATH%%R_LIBS%
Rcmd check plstools 
Rcmd build plstools
Rcmd INSTALL --build plstools
pause
