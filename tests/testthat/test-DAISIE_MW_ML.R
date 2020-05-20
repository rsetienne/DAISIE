context("DAISIE_MW_ML")

test_that("DAISIE_MW_ML produces correct output", {

  utils::data(archipelagos41)

M19_Nature_parameters<-c(4.007380e-02,  0.000000e+00,  1.945657e+00,  1.504297e-01,
                         Inf,  0.000000e+00,  6.725644e+01,
                          2.936351e-01,  5.909687e-02,  3.826885e-01,  2.651078e-02,
                         -3.653315e+03,  8.000000e+00,  0.000000e+00)

M19_computation<-DAISIE_MW_ML(
    datalist= archipelagos41,
    initparsopt = c(0.040073803,	1.945656546,	0.150429656,	67.25643672,	0.293635061,
                    0.059096872,	0.382688527, 0.026510781),
    idparsopt = c(1,3,4,7,8,9,10,11),
    parsfix = c(0,Inf,0) ,
    idparsfix = c(2,5,6),
    res = 100,
    ddmodel = 0,
    methode = 'ode45',
    cpus = 4,
    parallel = 'local',
    optimmethod = 'subplex',
    tol = c(1E-1, 1E-3, 1E-5),
    distance_type = 'continent',
    distance_dep = 'area_interactive_clado')


 expect_equal(M19_Nature_parameters,as.numeric(M19_computation),tol=0.000001)
})
