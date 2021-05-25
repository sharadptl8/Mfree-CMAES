function [obj,object]=objfn(pop,iter)
% iterations=1;
Parameter=10; % no. of parameters to optimize
lambda= 4+floor(3*log(Parameter)); % population size (i.e no. of sets)

% observed head values after 365 days
hobs=[60.11302004
57.0699172
54.05764992
53.1291622
53.16899608
66.18103521
63.4849499
57.97396408
55.71965028
55.14867121
54.4857801
53.25761192
54.29681475
54.41567478
76.15144758
74.32716905
69.88419729
63.1178414
58.82640086
58.90954283
58.85317534
57.95842239
56.88824794
85.24712509
81.73106157
77.59508807
69.18349594
65.62092545
65.87209647
65.07567268
63.21370081
61.41251031
57.40415756
90.23350472
89.5708635
86.52353327
80.16679679
77.36735342
74.87080035
72.17950529
68.89405628
66.2770624
63.98841687
63.24646008
100.1887102
99.05274883
96.88531579
92.30125412
87.95058152
83.31928114
79.17796119
74.70292899
70.72350079
67.19136843
110.3252829
107.8119259
104.9270532
101.125941
95.69715129
89.7757799
84.81051627
78.82063636
73.73047597
72.30479176
115.3627552
120.421298
117.7846721
115.2917256
111.294367
105.3711976
98.76482054
93.14620812
85.43384691
75.44364071
128.3441844
126.827347
124.7721797
120.5633655
114.935591
108.6304061
103.2042109
95.58035578
87.2466035
84.56491561
135.0551555
134.7441411
132.5771133
128.5181675
123.9132809
118.6972852
115.4104636
103.5344293
90.59193107
89.58743174
143.0141927
141.0331936
138.8140077
135.4561705
127.3862693
120.3584874
145.1513239
137.4356808
155.1139263
151.4034946
147.2967073
144.0987143
142.4391472
161.0785257
160.6940182
158.1909359
152.1181852
147.2387293
165.0785258
166.0785257
164.0785257
152.2451442
150.193327];

% noise=0+2.5*(randn(49,1)); % to generate random numbers of mean zero and
... std. deviation 2.5
    
% hobs=noise+hobs; % observed head values with noise

for p=1:lambda
    % Using the target vector (variable name "pop") from "DEalgo", head 
    ... values are simulated using FEM simulator by calling function call 
        ... "femdiff"
[hsimu]=Mfreediff(pop,p);

hsim(:,p)=hsimu(:,p); % saving the simulated head values "hsimu" with 
... different variable name "hsim"
    
for ob=1:117
    if  ob==1||ob==2||ob==8||ob==7||...
        ob==4||ob==12||ob==11||ob==33||...
        ob==24||ob==34||ob==65||ob==66||...ob==35||ob==36||...
        ob==29||ob==39||ob==58||ob==49||ob==60||...ob==50||ob==70||...
        ob==31||ob==41||ob==53||ob==63||...ob==64||ob==54||ob==74||...
        ob==80||ob==90||ob==100||ob==99||...ob==47||ob==68||...
        ob==93||ob==94||ob==92||ob==74||...ob==81||...
        ob==75||ob==85||ob==102||ob==89||...
        ob==97||ob==107||ob==112||ob==106||ob==111||...
        ob==101||ob==103||ob==104||ob==108||ob==109||ob==115
%         ob==106||ob==117||ob==103||ob==108||ob==110||ob==116%||ob==113||ob==114||ob==115||ob==104||ob==110||ob==106

    objvalue(ob,p)=((hsim(ob,p)-hobs(ob,:)).^2); % objective function 1 (squared error)
    % objvalue(ob,p)=abs(hsim(ob,p)-hobs(ob,:)); % objective function 2 (error)
    % objvalue(ob,p)=sqrt((((hsim(ob,p)-hobs(ob,:)).^2)/14)); % objective function 3 (root mean square error (RMSE))
    
end
end

obj(p,1)=sum(objvalue(:,p)); % add objective function values of all the observation wells

end
object(:,iter)=obj(:,1); % saves the objective function values of all the iterations correspondingly
end