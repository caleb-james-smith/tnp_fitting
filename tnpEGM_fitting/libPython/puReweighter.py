import ROOT as rt
import numpy as np
import sys
import argparse
import os

print '** puReweighter requires root_numpy.'
print '** To install on lxplus: '
print 'pip install --user root_numpy'
from root_numpy import  tree2array, array2tree


customWeights_17Nov2017MCv2 = {

    '2017_runB' : [ 2.78721e-09,
                    8.29537e-06,    0.0161464,      0.026547,       0.0470691,      0.0902818,      0.115647,       0.197857,
                    0.150448,       0.407495,       0.631471,       0.885343,       0.973195,       0.99756,        0.985337,
                    0.971295,       1.03354,        1.11933,        1.27721,        1.42536,        1.5472,         1.64752,
                    1.73905,        1.80084,        1.82923,        1.87493,        1.89668,        1.95769,        2.00439,
                    1.99327,        1.92761,        1.81131,        1.66118,        1.49815,        1.30579,        1.13063,
                    0.961473,       0.787239,       0.624429,       0.48327,        0.369995,       0.286518,       0.21152,
                    0.1527,         0.107546,       0.0723541,      0.0462112,      0.0272659,      0.0158567,      0.00874781,
                    0.00465276,     0.00245155,     0.00127945,     0.000665264,    0.000345037,    0.000174231,    9.2684e-05,
                    4.96701e-05,    2.68951e-05,    1.48571e-05,    8.49208e-06,    4.95017e-06,    2.93006e-06,    1.72019e-06,
                    1.05594e-06,    5.6375e-07,     3.00895e-07,    1.7454e-07,     1.06697e-07,    6.51644e-08,    3.67269e-08,
                    2.04583e-08,    8.83731e-09,    5.02539e-09,    1.47299e-09,    1.19777e-09,    3.39483e-10,    1.38009e-10,
                    3.51009e-11,    1.17901e-11,    1.48352e-12,    1.24697e-12,    4.94523e-13,    7.33467e-14,    2.21055e-14,
                    5.89623e-15,    4.46766e-15,    2.92693e-16,    1.06557e-16,    0,              0,              0,
                    0,              0,              0,              0,              0,              0,              0,
                    0,     ],
    '2017_runC' : [
         0.00026112,     0.0221736,      0.0128291,      0.0402351,      0.0831826,      0.0843967,      0.103655,
         0.162653,       0.110466,       0.219185,       0.260169,       0.437633,       0.713618,       1.07405,
         1.40912,        1.63428,        1.79933,        1.81333,        1.82617,        1.79701,        1.76691,
         1.75503,        1.76225,        1.75565,        1.72997,        1.73201,        1.72014,        1.75096,
         1.77686,        1.76028,        1.70267,        1.6046,         1.47877,        1.34319,        1.18282,
         1.03883,        0.899946,       0.753974,       0.614721,       0.491496,       0.391122,       0.317238,
         0.24766,        0.191283,       0.146115,       0.108253,       0.077375,       0.0519158,      0.0348403,
         0.0224469,      0.0140553,      0.00874558,     0.00537733,     0.00326879,     0.00195802,     0.00112472,
         0.000669624,    0.000395475,    0.000232925,    0.000138601,    8.48229e-05,    5.2811e-05,     3.34134e-05,
         2.10412e-05,    1.39333e-05,    8.08446e-06,    4.73116e-06,    3.03913e-06,    2.07966e-06,    1.43816e-06,
         9.28805e-07,    6.0019e-07,     3.04539e-07,    2.05993e-07,    7.27239e-08,    7.21124e-08,    2.52269e-08,
         1.28077e-08,    4.11458e-09,    1.76487e-09,    2.86572e-10,    3.13933e-10,    1.63882e-10,    3.2294e-11,
         1.30346e-11,    4.69725e-12,    4.78323e-12,    4.27695e-13,    2.12691e-13,    7.98385e-14,    1.24132e-13,
         2.83648e-14,    1.90417e-14,    9.3247e-15,     4.84836e-16,    1.78391e-16,    1.01492e-16,    0,
         0,              0,        ],


    '2017_runD' : [
        5.34841e-13,    1.06439e-05,    0.00754941,     0.0370866,      0.11429,        0.161627,       0.160145,
        0.193843,       0.0959119,      0.131594,       0.100569,       0.0861817,      0.0870812,      0.112852,
        0.202798,       0.456788,       0.952825,       1.41875,        1.71069,        1.85345,        1.93503,
        1.95719,        1.98401,        2.06412,        2.17753,        2.30768,        2.35005,        2.37909,
        2.34232,        2.19786,        1.96846,        1.68956,        1.40551,        1.14775,        0.907517,
        0.716872,       0.561438,       0.428375,       0.320232,       0.235508,       0.172047,       0.127177,
        0.0894371,      0.0613196,      0.0408963,      0.0259897,      0.0156401,      0.00866272,     0.00469865,
        0.00239154,     0.00115432,     0.000539369,    0.000242296,    0.000104587,    4.32024e-05,    1.66102e-05,
        6.42136e-06,    2.38797e-06,    8.58457e-07,    3.0211e-07,     1.05913e-07,    3.65733e-08,    1.24209e-08,
        4.06169e-09,    1.35058e-09,    3.80359e-10,    1.04387e-10,    3.03702e-11,    9.08708e-12,    2.65174e-12,
        6.97249e-13,    1.76923e-13,    3.39878e-14,    8.42428e-15,    1.00883e-15,    3.42283e-16,    1.134e-17,
        0,              0,              0,              0,              0,              0,              0,
        0,              0,              0,              0,              0,              0,              0,
        0,              0,              0,              0,              0,              0,              0,
        ],

    '2017_runE' : [
        9.976e-09,      0.00430579,     0.0197936,      0.0195893,      0.0469575,      0.0856007,      0.0898606,
        0.133994,       0.074991,       0.10824,        0.0873948,      0.118887,       0.178507,       0.249594,
        0.318471,       0.389619,       0.490148,       0.593916,       0.727234,       0.848012,       0.941199,
        1.00838,        1.05839,        1.07527,        1.05638,        1.03848,        1.0103,         1.01924,
        1.04871,        1.08311,        1.12128,        1.15404,        1.17759,        1.19714,        1.19476,
        1.21012,        1.23627,        1.25266,        1.26754,        1.28955,        1.33681,        1.44395,
        1.5323,         1.6388,         1.76053,        1.85531,        1.89803,        1.82405,        1.74553,
        1.58955,        1.39009,        1.19177,        0.995484,       0.810875,       0.642734,       0.483302,
        0.373428,       0.284424,       0.215222,       0.164329,       0.129213,       0.103762,       0.0852167,
        0.070273,       0.0616186,      0.0479701,      0.0382311,      0.0339894,      0.0327391,      0.0324185,
        0.0304894,      0.0291627,      0.0222454,      0.022956,       0.0125378,      0.0194903,      0.0108263,
        0.0088362,      0.00461907,     0.00326246,     0.000882622,    0.00162985,     0.00145087,     0.000493156,
        0.000347263,    0.000220807,    0.000401183,    6.47169e-05,    5.87115e-05,    4.06217e-05,    0.000117842,
        5.07655e-05,    6.49354e-05,    6.13054e-05,    6.33085e-06,    4.22432e-06,    6.19473e-06,    3.46389e-06,
        5.66718e-06,    1.1995e-06,
        ],
    '2017_runF' : [
         0.000791184,    0.113765,       0.161239,       0.134982,       0.09516,        0.130148,       0.116681,
         0.115613,       0.131078,       0.753048,       1.01439,        1.24964,        1.06256,        0.783467,
         0.58507,        0.488352,       0.475015,       0.465189,       0.474411,       0.48612,        0.525865,
         0.604733,       0.694277,       0.743258,       0.742267,       0.730203,       0.703604,       0.694884,
         0.693959,       0.69485,        0.703141,       0.718069,       0.738532,       0.767131,       0.791625,
         0.835908,       0.891396,       0.937885,       0.980591,       1.03486,        1.13255,        1.32895,
         1.58039,        1.93993,        2.4177,         2.94552,        3.43089,        3.66776,        3.79819,
         3.63427,        3.24338,        2.76013,        2.23112,        1.72014,        1.26752,        0.874813,
         0.616495,       0.428977,       0.299926,       0.216405,       0.166374,       0.136477,       0.120221,
         0.111475,       0.114341,       0.107111,       0.10444,        0.114357,       0.135551,       0.164229,
         0.18735,        0.215157,       0.194953,       0.236431,       0.150195,       0.268901,       0.170417,
         0.157279,       0.092175,       0.0723903,      0.0216036,      0.0436668,      0.0422279,      0.0154774,
         0.0116665,      0.0078836,      0.0151136,      0.00255425,     0.00241046,     0.00172259,     0.00512494,
         0.00224817,     0.0029075,      0.00275556,     0.000283618,    0.000187269,    0.000269796,    0.000147143,
         0.000233107,    4.74285e-05,
         ],

    '2017_runBCDEF' : [0.000319471, 0.043333,  0.0627865,  0.0647238,  0.0780242, 0.108199,  0.111983,
                 0.148095,    0.112399,  0.381801,   0.494781,   0.647735,  0.674107,  0.687767,
                 0.723788,    0.784618,  0.898843,   0.979519,   1.0633,    1.11899,   1.16819,
                 1.22006,     1.27549,   1.30905,    1.31339,    1.32447,   1.31354,   1.32983,
                 1.34375,     1.33185,   1.29874,    1.24632,    1.18278,   1.11995,   1.04354,
                 0.988411,    0.944833,  0.896181,   0.851375,   0.820549,  0.820257,  0.877121,
                 0.9505,      1.06921,   1.23497,    1.41449,    1.57161,   1.62378,   1.64322,
                 1.5508,      1.37578,   1.1718,     0.953895,   0.744876,  0.558862,  0.394524,
                 0.285334,    0.20402,   0.146289,   0.107565,   0.0833341, 0.0678563, 0.0583878,
                 0.0521837,   0.0511854, 0.0457683,  0.042717,   0.0450134, 0.0516682, 0.0609794,
                 0.0681131,   0.0769145, 0.068758,   0.082488,   0.0519442, 0.0923373, 0.0581781,
                 0.0534347,   0.031191,  0.0244144,  0.00726572, 0.0146515, 0.0141407, 0.00517426,
                 0.00389476,  0.00262877,0.00503468, 0.000850185,0.00080179,0.00057269,0.00170316,
                 0.000746922, 0.00096581,0.00091529, 9.42116e-05,6.2216e-05,8.9658e-05,4.89178e-05,
                 7.75365e-05, 1.57862e-05],

    }

puMC = {
    'Spring2016MC_PUscenarioV1' : [ 0.000829312873542, 0.00124276120498, 0.00339329181587, 0.00408224735376, 0.00383036590008, 
                                    0.00659159288946,  0.00816022734493, 0.00943640833116, 0.0137777376066,  0.017059392038,
                                    0.0213193035468,   0.0247343174676,  0.0280848773878,  0.0323308476564,  0.0370394341409,  
                                    0.0456917721191,   0.0558762890594,  0.0576956187107,  0.0625325287017,  0.0591603758776,
                                    0.0656650815128,   0.0678329011676,  0.0625142146389,  0.0548068448797,  0.0503893295063,  
                                    0.040209818868,    0.0374446988111,  0.0299661572042,  0.0272024759921,  0.0219328403791,
                                    0.0179586571619,   0.0142926728247,  0.00839941654725, 0.00522366397213, 0.00224457976761, 
                                    0.000779274977993, 0.000197066585944,7.16031761328e-05,0.0             , 0.0,
                                    0.0,        0.0,        0.0,        0.0,        0.0,    
                                    0.0,        0.0,        0.0,        0.0,        0.0],
    
    'Moriond17MC_mix_2016'      : [ 1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,
                                    0.00071209  , 0.00130121 ,0.00245255  ,0.00502589  ,0.00919534  ,0.0146697   ,0.0204126   ,
                                    0.0267586   ,0.0337697   ,0.0401478   ,0.0450159   ,0.0490577   ,0.0524855   ,0.0548159   ,
                                    0.0559937   ,0.0554468   ,0.0537687   ,0.0512055   ,0.0476713   ,0.0435312   ,0.0393107   ,
                                    0.0349812   ,0.0307413   ,0.0272425   ,0.0237115   ,0.0208329   ,0.0182459   ,0.0160712   ,
                                    0.0142498   ,0.012804    ,0.011571    ,0.010547    ,0.00959489  ,0.00891718  ,0.00829292  , 
                                    0.0076195   ,0.0069806   ,0.0062025   ,0.00546581  ,0.00484127  ,0.00407168  ,0.00337681  ,
                                    0.00269893  ,0.00212473  ,0.00160208  ,0.00117884  ,0.000859662 ,0.000569085 ,0.000365431 ,
                                    0.000243565 ,0.00015688  ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05  ,
                                    1.73032e-05 ,1.435e-05   ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,
                                    1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05 ],

   
    'mix_2017_25ns_UltraLegacy_PoissonOOTPU' : [ 1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
                                                 0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
                                                 0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
                                                 0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
                                                 0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
                                                 0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
                                                 0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
                                                 0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
                                                 0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
                                                 0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
                                                 0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
                                                 0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
                                                 0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
                                                 3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
                                                 1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0],


    'Moriond18MC_mix_2017' :[3.39597497605e-05,    6.63688402133e-06,     1.39533611284e-05,     3.64963078209e-05,     6.00872171664e-05,     9.33932578027e-05,    0.000120591524486,
                             0.000128694546198,    0.000361697233219,     0.000361796847553,     0.000702474896113,     0.00133766053707,      0.00237817050805,     0.00389825605651,
                             0.00594546732588,     0.00856825906255,      0.0116627396044,       0.0148793350787,       0.0179897368379,       0.0208723871946,      0.0232564170641,
                             0.0249826433945,      0.0262245860346,       0.0272704617569,       0.0283301107549,       0.0294006137386,       0.0303026836965,      0.0309692426278,   
                             0.0308818046328,      0.0310566806228,       0.0309692426278,       0.0310566806228,       0.0310566806228,       0.0310566806228,      0.0307696426944,
                             0.0300103336052,      0.0288355370103,       0.0273233309106,       0.0264343533951,       0.0255453758796,       0.0235877272306,      0.0215627588047,
                             0.0195825559393,      0.0177296309658,       0.0160560731931,       0.0146022004183,       0.0134080690078,       0.0129586991411,      0.0125093292745,
                             0.0124360740539,      0.0123547104433,       0.0123953922486,       0.0124360740539,       0.0124360740539,       0.0123547104433,      0.0124360740539,
                             0.0123387597772,      0.0122414455005,       0.011705203844,        0.0108187105305,       0.00963985508986,      0.00827210065136,     0.00683770076341,
                             0.00545237697118,     0.00420456901556,      0.00367513566191,      0.00314570230825,      0.0022917978982,       0.00163221454973,     0.00114065309494,
                             0.000784838366118,    0.000533204105387,     0.000358474034915,     0.000238881117601,     0.0001984254989,       0.000157969880198,    0.00010375646169,
                             6.77366175538e-05,    4.39850477645e-05,     2.84298066026e-05,     1.83041729561e-05,     1.17473542058e-05,     7.51982735129e-06,    6.16160108867e-06,
                             4.80337482605e-06,    3.06235473369e-06,     1.94863396999e-06,     1.23726800704e-06,     7.83538083774e-07,     4.94602064224e-07,    3.10989480331e-07,
                             1.94628487765e-07,    1.57888581037e-07,     1.2114867431e-07,      7.49518929908e-08,     4.6060444984e-08,      2.81008884326e-08,    1.70121486128e-08,
                             1.02159894812e-08],


#https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi.py

    'mix_2018_25ns_UltraLegacy_PoissonOOTPU' :[8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
   0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
    0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
    0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
    0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
    0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
    0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
    0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
    0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
    0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
    0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
    0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
    0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
    2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
    3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
    4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0],

    ## copied from https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_UltraLegacy_PoissonOOTPU_cfi.py
    'mix_2016_25ns_UltraLegacy_PoissonOOTPU' :[1.00402360149e-05, 5.76498797172e-05, 7.37891400294e-05, 0.000110932895295, 0.000158857714773,
   0.000368637432599, 0.000893114107873, 0.00189700774575, 0.00358880167437, 0.00636052573486,
    0.0104173961179, 0.0158122597405, 0.0223785660712, 0.0299186888073, 0.0380275944896,
    0.0454313901624, 0.0511181088317, 0.0547434577348, 0.0567906239028, 0.0577145461461,
    0.0578176902735, 0.0571251566494, 0.0555456541498, 0.053134383488, 0.0501519041462,
    0.0466815838899, 0.0429244592524, 0.0389566776898, 0.0348507152776, 0.0307356862528,
    0.0267712092206, 0.0229720184534, 0.0193388653099, 0.0159602510813, 0.0129310510552,
    0.0102888654183, 0.00798782770975, 0.00606651703058, 0.00447820948367, 0.00321589786478,
    0.0022450422045, 0.00151447388514, 0.000981183695515, 0.000609670479759, 0.000362193408119,
    0.000211572646801, 0.000119152364744, 6.49133515399e-05, 3.57795801581e-05, 1.99043569043e-05,
    1.13639319832e-05, 6.49624103579e-06, 3.96626216416e-06, 2.37910222874e-06, 1.50997403362e-06,
    1.09816650247e-06, 7.31298519122e-07, 6.10398791529e-07, 3.74845774388e-07, 2.65177281359e-07,
    2.01923536742e-07, 1.39347583555e-07, 8.32600052913e-08, 6.04932421298e-08, 6.52536630583e-08,
    5.90574603808e-08, 2.29162474068e-08, 1.97294602668e-08, 1.7731096903e-08, 3.57547932012e-09,
    1.35039815662e-09, 8.50071242076e-09, 5.0279187473e-09, 4.93736669066e-10, 8.13919708923e-10,
    5.62778926097e-09, 5.15140589469e-10, 8.21676746568e-10, 0.0, 1.49166873577e-09,
    8.43517992503e-09, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0],
}

### MC pu scenario to be used
#puMCscenario = 'Spring2016MC_PUscenarioV1'
#puMCscenario = 'mix_2017_25ns_UltraLegacy_PoissonOOTPU'
#puMCscenario = 'mix_2018_25ns_UltraLegacy_PoissonOOTPU'
puMCscenario = 'mix_2016_25ns_UltraLegacy_PoissonOOTPU'
customWeightsName= 'weights'
###puDirEOS = '/eos/cms/store/group/phys_egamma/swmukher/ntuple_2017/PU/'
#puDirEOS = '/eos/cms/store/group/phys_egamma/soffi/TnP/ntuples_01162018/PU/'
#puDirEOS = '/eos/cms/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2018_PUData/'
puDirEOS = '/eos/cms/store/group/phys_egamma/asroy/Tag-and-Probe_Tree/UL2016_PUData/'


#### Compute weights for all data epoch specified below
puDataEpoch = {
    '2016_runBCD' : puDirEOS + 'pu_dist_runBCD_692.root',
    '2016_runEF' : puDirEOS + 'pu_dist_runEF_692.root',
    '2016_runGH' : puDirEOS + 'pu_dist_runGH_692.root',
    '2016_run2016' : puDirEOS + 'pu_dist_run2016_692.root',
    }
'''
puDataEpoch = {
    '2018_runA' : puDirEOS + 'pileup_2018_RunA.root',
    '2018_runB' : puDirEOS + 'pileup_2018_RunB.root',
    '2018_runC' : puDirEOS + 'pileup_2018_RunC.root',
    '2018_runD'  : puDirEOS +'pileup_2018_RunD.root' ,
    '2018_runABCD' : puDirEOS + 'pileup_2018_RunABCD.root',
    }
puDataEpoch = {
    '2017_runB' : puDirEOS + 'pileup_2017_RUNB.root',
    '2017_runC' : puDirEOS + 'pileup_2017_RUNC.root',
    '2017_runD'  : puDirEOS +'pileup_2017_RUND.root' ,
    '2017_runE'  : puDirEOS +'pileup_2017_RUNE.root' ,
    '2017_runF' : puDirEOS + 'pileup_2017_RUNF.root',    
    '2017_runBCDEF' : puDirEOS + 'pileup_2017_41fb.root',
    }
'''
nVtxDataEpoch = {
    '2016_runBCD' : 'etc/inputs/nVtx_2016_runBCD.root',
    '2016_runEF'  : 'etc/inputs/nVtx_2016_runEF.root' ,
    '2016_runGH'  : 'etc/inputs/nVtx_2016_runGH.root' ,
}

rhoDataEpoch = {
    '2016_runE'   : 'etc/inputs/rho_pu_runE.root',
    '2016_runGH'  : 'etc/inputs/rho_pu_runGH.root',
}





def reweight( sample, puType = 0,useCustomW=False  ):
    if sample.path is None:
        print '[puReweighter]: Need to know the MC tree (option --mcTree or sample.path)'
        sys.exit(1)
    

### create a tree with only weights that will be used as friend tree for reweighting different lumi periods
    print 'Opening mc file: ', sample.path[0]
    fmc = rt.TFile(sample.path[0],'read')
    tmc = None
    if sample.tnpTree is None:
        dirs = fmc.GetListOfKeys()
        for d in dirs:
            if (d.GetName() == "sampleInfo"): continue
            tmc = fmc.Get("%s/fitter_tree" % d.GetName())
    else:
        tmc = fmc.Get(sample.tnpTree)
    

#### can reweight vs nVtx but better to reweight v truePU
    puMCnVtx = []
    puMCrho = []
    if   puType == 1 :
        hmc   = rt.TH1F('hMC_nPV'  ,'MC nPV'  , 75,-0.5,74.5)
        tmc.Draw('event_nPV>>hMC_nPV','','goff')
        hmc.Scale(1/hmc.Integral())
        for ib in range(1,hmc.GetNbinsX()+1):
            puMCnVtx.append( hmc.GetBinContent(ib) )
        print 'len nvtxMC = ',len(puMCnVtx)

    elif puType == 2 :
        hmc   = rt.TH1F('hMC_rho'  ,'MC #rho'  , 75,-0.5,74.5)
        tmc.Draw('rho>>hMC_rho','','goff')
        hmc.Scale(1/hmc.Integral())
        for ib in range(1,hmc.GetNbinsX()+1):
            puMCrho.append( hmc.GetBinContent(ib) )
        print 'len rhoMC = ',len(puMCrho)
    

    puDataDist = {}
    puDataArray= {}
    weights = {}
    epochKeys = puDataEpoch.keys()
    if puType == 1  : epochKeys = nVtxDataEpoch.keys()
    if puType == 2  : epochKeys = rhoDataEpoch.keys()
 
    for pu in epochKeys:
        fpu = None
        if   puType == 1 : fpu = rt.TFile(nVtxDataEpoch[pu],'read')
        elif puType == 2 : fpu = rt.TFile(rhoDataEpoch[pu],'read')
        else             : fpu = rt.TFile(puDataEpoch[pu],'read')
        puDataDist[pu] = fpu.Get('pileup').Clone('puHist_%s' % pu)
        puDataDist[pu].Scale(1./puDataDist[pu].Integral())
        puDataDist[pu].SetDirectory(0)
        puDataArray[pu] = []
        for ipu in range(len(puMC[puMCscenario])):
            ibin_pu  = puDataDist[pu].GetXaxis().FindBin(ipu+0.00001)
            puDataArray[pu].append(puDataDist[pu].GetBinContent(ibin_pu))
        print 'puData[%s] length = %d' % (pu,len(puDataArray[pu]))
        fpu.Close()
        weights[pu] = []

    mcEvts = tree2array( tmc, branches = ['weight','truePU','event_nPV','rho'] )


    pumc = puMC[puMCscenario]
    if   puType == 1:  pumc = puMCnVtx
    elif puType == 2:  pumc = puMCrho
    else            :  pumc = puMC[puMCscenario]

    puMax = len(pumc)
    print '-> nEvtsTot ', len(mcEvts)
#    print "--------------------------" 
    for ievt in xrange(len(mcEvts)):
        if ievt%100000 == 0 :            print 'iEvt:',ievt
#        print 'iEvt:',ievt
        evt = mcEvts[ievt]
        for pu in epochKeys:
#            print pu
#AR            customWeights=customWeights_17Nov2017MCv2[pu]
#            print customWeights_17Nov2017MCv2[pu]
#            print "--------------------------"
            pum = -1
            pud = -1
            if puType == 1 and evt['event_nPV'] < puMax:
                pud = puDataArray[pu][evt['event_nPV']]
                pum = pumc[evt['event_nPV']]
            if puType == 2 and int(evt['rho']) < puMax:
                pud = puDataArray[pu][int(evt['rho'])]
                pum = pumc[int(evt['rho'])]
            elif puType == 0:
#                if ievt%1000: print pu, evt['truePU']
                if evt['truePU']> 0 and evt['truePU']<99: 
                    pud = puDataArray[pu][evt['truePU']] 
                    pum = pumc[evt['truePU']]
            puw = 0
            if pum > 0: 
                puw  = pud/pum
#                if use customized weights
            if useCustomW:
                puw=0
                if  evt['truePU']> 0 and evt['truePU']<97:
                    puw = customWeights[evt['truePU']]
                
#                    print evt['truePU'],puw

            if evt['weight'] > 0 : totw = +puw
            else                 : totw = -puw
            weights[pu].append( ( puw,totw) )
#    print "====================="
#    print weights[pu]

    newFile    = rt.TFile( sample.puTree, 'recreate')

    for pu in epochKeys:
        treeWeight = rt.TTree('weights_%s'%pu,'tree with weights')
        wpuarray = np.array(weights[pu],dtype=[('PUweight',float),('totWeight',float)])
        array2tree( wpuarray, tree = treeWeight )
        treeWeight.Write()

    newFile.Close()    
    fmc.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='tnp EGM pu reweighter')
    parser.add_argument('--mcTree'  , dest = 'path',  default = None, help = 'MC tree to compute weights for')
    parser.add_argument('puTree'    , default = None                , help = 'output puTree')

    args = parser.parse_args()
    args.path = [args.path]
    reweight(args)




