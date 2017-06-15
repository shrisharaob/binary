(TeX-add-style-hook
 "p7g0_m75_T10"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("inputenc" "utf8") ("geometry" "left=1cm" "right=1cm" "top=.8cm" "bottom=1.2cm")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "fontenc"
    "inputenc"
    "lmodern"
    "geometry"
    "graphicx"
    "subcaption")))

