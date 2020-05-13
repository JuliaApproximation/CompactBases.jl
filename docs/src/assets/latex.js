requirejs.config({
    paths: {
        'mathjax': 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS_HTML',
    },
    shim: {
        'mathjax' : {
            exports: "MathJax"
        },
    }
});

require(['mathjax'], function(MathJax) {
    MathJax.Hub.Config({
        TeX: {
            Macros: {
                defd: "≝",
                ket: ["|#1\\rangle",1],
                bra: ["\\langle#1|",1],
                braket: ["\\langle#1|#2\\rangle",2],
                ketbra: ["|#1\\rangle\\!\\langle#2|",2],
                matrixel: ["\\langle#1|#2|#3\\rangle",3],
                vec: ["\\mathbf{#1}",1],
                mat: ["\\mathsf{#1}",1],
                conj: ["#1^*",1],
                im: "\\mathrm{i}",
                operator : ["\\mathfrak{#1}",1],
                Hamiltonian : "\\operator{H}",
                hamiltonian : "\\operator{h}",
                Lagrangian : "\\operator{L}",
                fock : "\\operator{f}",
                lagrange : ["\\epsilon_{#1}",1],
                vary : ["\\delta_{#1}",1],
                onebody : ["(#1|#2)",2],
                twobody : ["[#1|#2]",2],
                twobodydx : ["[#1||#2]",2],
                direct : ["{\\operator{J}_{#1}}",1],
                exchange : ["{\\operator{K}_{#1}}",1],
                diff : ["\\mathrm{d}#1\\,",1],
                B : ["\\mathrm{B}_{#1,#2}", 2],
                bmat : ["\\begin{bmatrix}#1\\end{bmatrix}", 1],
                ceil : ["\\left\\lceil #1\\right\\rceil", 1],
                floor : ["\\left\\lfloor #1\\right\\rfloor", 1],
                space : ["\\mathcal{#1}", 1]
            }
        }
    });
})
