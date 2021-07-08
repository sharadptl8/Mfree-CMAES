function [hsimu]= Mfreediff(pop,p)
nnodes=117;         % No of nodes throughout the domian
alpha=3;            % Spport to the RBF               
dc=5389.32;         % Avg. Nodal spacing
Cs=alpha*dc;        % Shape parameter
Aera= 2798950000;   % Area of MRBC Aquifer
dt=1;               % Time step size 

% X-distance is specified as column vector (i.e array of size 117 x 1)

% Intial head values for each node
h=[60 
57
54
53
53
66
64
50
57
54
53
53
54
54
76
75
67
63
50
58
56
55
54
85
80
75
65
63
63
63
60
60
57
90
87
83
75
73
70
68
66
64
63
63
100
95
93
86
84
79
75
73
70
67
110
105
100
95
90
86
83
76
74
72
115
120
114
110
107
102
95
90
83
75
128
125
120
115
110
105
100
93
85
84
135
133
130
125
123
118
115
103
90
88
143
140
138
135
127
120
145
137
155
150
147
144
142
161
160
158
150
147
165
166
164
152
150];

for l=1:117   
K(117,1)=0;% Nodal hydraukic conductivity matrix in X-direction
if l==1||l==2||l==3||l==6||l==7||l==8||l==9||l==10||l==15||l==16||l==17||l==18||l==19|| l==20
    K(l,:)=K(l,:)+pop(1,p); % extracting first row from population 
elseif l==101||l==103||l==104||l==108||l==109||l==110||l==113||l==114||l==115||l==116||l==117
    K(l,:)=K(l,:)+pop(10,p);
end
end
%% Initialization of diffirent matrix in coefficient matrix 

R=zeros(117,117); % RBF matrix
KD2=zeros(117,117); % Second order derivative of RBF wrt to y
G=zeros(117,117);% Coefficient matrix for unsteady state

for j=1:117
     for nnode=1:117
        R(nnode,j)=(sqrt(abs((x(nnode)-x(j))^2+(y(nnode)-y(j))^2)+Cs^2));
    end
end

T=bsxfun(@times,R,Sy);
hKx=h.*K;
hKy=h.*K;

dt=1;  
for i=1:117 
        for j=1:117
          if i==1||i==2||i==3||i==4||i==5||i==12||i==13||i==14||i==33||i==44||i==54||i==64||i==74||...
            i==84||i==94||i==93||i==92||i==91||i==100||i==99||i==102||i==107||i==112||i==117||i==116||i==110||... 
             i==115||i==114||i==113||i==108||i==103||i==101||i==95||i==85||i==75||i==66||i==65||i==55||i==45||... 
             i==34||i==24||i==15||i==6
            G(i,j)=sqrt(abs((x(i)-x(j))^2+(y(i)-y(j))^2)+Cs^2);% boundary condition will be same as in steady state 
          else 
             G(i,j)=(T(i,j)/dt)-((((y(i)-y(j))^2)+(Cs^2))/((sqrt(abs((x(i)-x(j))^2+(y(i)-y(j))^2)+Cs^2))^3)*hKx(i))-((((x(i)-x(j))^2)+(Cs^2))/((sqrt(abs((x(i)-x(j))^2+(y(i)-y(j))^2)+Cs^2))^3)*hKy(i));
          end
        end        
end  
for z=1:365
for M=1:117
for l=1:117
    if  M==1||M==2||M==3||M==4||M==5||M==12||M==13||M==14||M==33||M==44||M==54||M==64||M==74||...
        M==84||M==94||M==93||M==92||M==91||M==100||M==99||M==102||M==107||M==112||M==117||M==116||M==110||... 
        M==115||M==114||M==113||M==108||M==103||M==101||M==95||M==85||M==75||M==66||M==65||M==55||M==45||... 
        M==34||M==24||M==15||M==6
        Fadd(M,1)=h(M,1);
    else 
        Fadd(M,1)=(((L(M,1)).*(Sy(M,1)))/(R(M,l)))/dt;
    end
        Ffinal(M,1)=Su(M,1)+Fadd(M,1); 
end
end  
 h=R*(G\Ffinal);
 head=h;
 hsimu(:,p)=head;
end 
end
