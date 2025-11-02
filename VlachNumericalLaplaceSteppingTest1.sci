/*
   TIME DOMAIN RESPONSE OF TRANSFER FUNCTIONS (TF) BY NUMERICAL
   INVERSION OF THE LAPLACE TRANSFORM.
   THIS PROGRAM WILL COMPUTE THE IMPULSE RESPONSE OF A TRANSFER
   FUNCTION GIVEN IN THE FORM:
             B(1) + B(2)*S + ... + B(M)*S**(M-1)
          ------------------------------------------
          A(1) + A(2)*S + ,.. + A(N)*S**(N-1) + S**N

   WITH N.GE.M , N,LE.15 , AND N.GE.2 IN THIS IMPLEMENTATION
   THE SOLUTION IS FOUND TILL TIME TEND WITH AN INTERNAL STEP OF
   DELT AND WITH OUTPUT PRODUCED EVERY outKK*DELT TIME UNITS

   M,N selection: p=M+N ; p=highest order of denominator.
   Select N = M-2 for all practical cases.
   For all practical purposes stay close to M+N=p !

   Theory in:
     Method for Computing Time Response of Systems Described
     by Transfer Function
     K.Singhal, J.Vlach
     Journal of the Franklin Institute, 1981

   Please be advised that the program uses a unit step pulse as
   input function per default.
   If you want to tinker with it, try it
   below @ "input function W(s)" in the Time stepping loop.

   Heiko Schroeter, University of Bremen, Nov 2025
*/
clc; // clear console
clf; // clear plot window
banner();
clearglobal();

// predefine size of denominator and nominator polynomial
A = zeros(16,1);
B = zeros(15,1);

// Transfer Function setup:
// In the polynom A (denominator) the
// term s^N is "build in" and therefore
// not set in Array A below.
M      = 2;
N      = 2;
B(1:M) = [25,1];    // B(1)+B(2)*s+...+B(M)*s**(M-1)
A(1:N) = [5,0.5];   // A(1)+A(2)*s+...+A(N)*s**(N-1) + "s**N"
// If you want at least having a few more Zeros than power of
// denominator, add "zeroPrec=NNN" to N
// for Zero/Residue selection.
// Check results for accuracy !
zeroPrec = 0;

// Time:
h      =  0.2; // Stepsize for numerical TF
TEND   = 20.0; // End time
KK     =  0.2; // Output every t*KK times

// Collect zeros and residues in a list to avoid excesive
// powers of Zeros and Residues during calculations.
// We select proper Zeros/Residues from the
// power of N in this program, hopefully ...
nullst=list();

// Pade Zero/Residue library.
// M=24,N=22
Z24=[complex(30.4053703404472890,1.6966885646518453),
     complex(30.1733650328535300,5.0060695539676909),
     complex(29.6623905107852600,8.4071674777829060),
     complex(28.8486156353132740,11.7732451193476600),
     complex(27.7770463258912680,15.1696779756105850),
     complex(26.3787445913999870,18.6041674001305860),
     complex(24.6325937075632860,22.0770631023638000),
     complex(22.4843550296975440,25.6133397353179900),
     complex(19.8502490392111400,29.2401445373629500),
     complex(16.5924461730515700,33.0071828528970900),
     complex(12.4455857549677450,37.0153908423283970),
     complex(6.7317632551660800,41.5429572419370000)];
KP24=[complex(2433391368255.78680,-17418545581916.062),
      complex(-5831324107448.0439,12433328721581.7190),
      complex(5496816340442.18000,-6068129779319.8799),
      complex(-3292806642424.6230,1923543951868.90050),
      complex(1313153797125.06060,-248441383293.35876),
      complex(-338738057374.48083,-61441452839.569763),
      complex(52755856940.5937200,35394751214.9989200),
      complex(-3838822855.4605732,-6960970355.9977188),
      complex(-40160272.367098898,623584264.758136000),
      complex(19022592.5201458900,-19377561.121864825),
      complex(-505460.74786890252,-83608.486733299855),
      complex(-191.54123066266891,2303.51715239361160)];
// M=22,N=20
Z22 =[complex(27.7758272188982160,1.6708876804019908),
      complex(27.4947799464109060,5.0160503975764378),
      complex(26.9264909552433500,8.3706317381326090),
      complex(26.0579121569180870,11.7423272687191290),
      complex(24.8672085299304000,15.1401818966264960),
      complex(23.3210860094199450,18.5763065284016960),
      complex(21.3687102868970000,22.0681351477798680),
      complex(18.9302992287510800,25.6427848047375180),
      complex(15.8714505516824930,29.3473310274765300),
      complex(11.9332029149560290,33.2776968620175400),
      complex(6.4530671439250970,37.6994854291197000)];
KP22=[complex(174223125639.053660,-1218345158935.5286),
      complex(-390207246922.42212,852472842262.469600),
      complex(359232906090.955400,-404344126157.92200),
      complex(-201376087165.63409,117930391776.165620),
      complex(72473317437.2514700,-13813457561.641846),
      complex(-16255805808.024391,-3166723764.2094526),
      complex(2033046076.82846200,1483198426.50774960),
      complex(-97475816.538217857,-219644196.43236592),
      complex(-3007405.6302386597,12376350.4240420020),
      complex(267395.573235217900,-127463.52255101516),
      complex(-988.60059864539085,-1372.4680715780996)];
// M=20,N=18
Z20 =[complex(25.1216083683247980,1.6646223827928038),
      complex(24.8126465458755400,4.9973217333814630),
      complex(24.1873298400001650,8.3413747653076560),
      complex(23.2272965499030070,11.7063106597156340),
      complex(21.9027613447091400,15.1034139501942950),
      complex(20.1663314802179850,18.5493539972803040),
      complex(17.9416571508178440,22.0696839499936500),
      complex(15.0987329137105640,25.7086962305363770),
      complex(11.3855752704304760,29.5569516096536490),
      complex(6.1561513477865720,33.8664317222223000)];
KP20=[complex(12021511652.2597320,-85411443711.57316600),
      complex(-26302647433.16952500,58389784585.643770000),
      complex(23059116318.537626000,-26365599892.64752200),
      complex(-11933819591.32598900,7083953039.8706110000),
      complex(3806800685.7747910000,-716286610.4397175300),
      complex(-711865935.8487265100,-154712648.5875087400),
      complex(66290275.359154180000,54767678.286548059000),
      complex(-1571414.589567282500,-5429194.513035957700),
      complex(-83548.02339305639900,141956.63072781426000),
      complex(1173.7378344646889000,298.4019627022859900)];
// M=18,N=16
Z18 =[complex(22.4664015786950600,1.6565077400455628),
      complex(22.1240565688110040,4.9745359953432120),
      complex(21.4281601739241640,8.3060524059258690),
      complex(20.3541531623414480,11.6632860713996470),
      complex(18.8587732727316440,15.0616707719535730),
      complex(16.8704180971027000,18.5251536010159170),
      complex(14.2642504525815820,22.0951208604098470),
      complex(10.7961224309839390,25.8562752851763300),
      complex(5.8377624346713810,30.0457311819202200)];
KP18=[complex(825575812.4087818000,-5979378520.58621880),
      complex(-1755752721.26384070,3973391761.0805799),
      complex(1450261500.749738800,-1689259761.6161900),
      complex(-680261576.619928840,409298093.105535900),
      complex(186211906.4119433600,-33782668.1142548400),
      complex(-27342631.2837899630,-6980987.33409900590),
      complex(1653435.386945674500,1666411.19815059500000),
      complex(-3925.20500703081600,-87347.5725929117180000),
      complex(-765.445714966286350,357.9153018508356700)];
// M=16,N=14
Z16 =[complex(19.8101574592180660,1.6468145348069050),
      complex(19.4260023882312440,4.9459053647182948),
      complex(18.6414951183484890,8.2623328097197000),
      complex(17.4205593558612100,11.6112396645079770),
      complex(15.6976183438144460,15.0149789161144750),
      complex(13.3544760009097080,18.5117397680788340),
      complex(10.1560781516918750,22.1798376257739790),
      complex(5.4936131860387360,26.2399655049469880)];
KP16=[complex(56256592.65559971000,-417742932.32983762),
      complex(-115635597.169700240,268062690.69494219),
      complex(88726603.37493292000,-105629202.14855157),
      complex(-36732668.3626251740,22372243.916864687),
      complex(8215579.521955218000,-1368661.5003691358),
      complex(-855586.56994330627,-277226.20656531898),
      complex(24822.7487324780280,36168.85760441191),
      complex(253.962585329660760,-511.86018939523564)];
// M=14,N=12
Z14 =[complex(17.1526203036969000,1.6341762131881524),
      complex(16.7147257848737050,4.9093779970759870),
      complex(15.8146338580616120,8.2067188334503140),
      complex(14.3962858335385320,11.5471965661553150),
      complex(12.3505043646558320,14.9656842613605200),
      complex(9.4533579165358260,18.5334197604499100),
      complex(5.1178718946131310,22.4527018231580240)];
KP14=[complex(3794630.02454944350,-29101101.2235581690),
      complex(-7471287.3435699772,17862919.17846260400),
      complex(5219303.66986495390,-6378881.17030751980),
      complex(-1832796.0566160935,1124727.107025519300),
      complex(308664.009573355980,-41731.0263554563860),
      complex(-18587.352203621445,-8687.80748878718800),
      complex(72.2886029684291500,364.5905114500656000)];
// M=12,N=10
Z12 =[complex(14.4929942543921480,1.6173310060964020),
      complex(13.9836392308959460,4.8608110270581900),
      complex(12.9259172792931630,8.1336223787795080),
      complex(11.2246068338888730,11.4673703426005520),
      complex(8.6705609733887200,14.9254273078465990),
      complex(4.7022814281112599,18.6890942420081800)];
KP12=[complex(251895.719679733860,-2018607.55517570230),
      complex(-469237.91613057588,1168946.41052460200),
      complex(289569.669958283200,-365339.932857681940),
      complex(-80809.075539924990,49317.39144909328000),
      complex(8751.56856924547600,-566.4028686922704300),
      complex(-169.96656637631784,-155.0735720910517400)];
// M=10,N=8
Z10 =[complex(11.83009373916912  , 1.593753005885994),
      complex(11.22085377939504  , 4.7929641675657),
      complex( 9.933383722174419 , 8.033106334266357),
      complex( 7.781146264465029 ,11.36889164904975),
      complex( 4.23452249479694  ,14.95704378128163)];
KP10=[complex( 16286.62368059556, -139074.71155151649),
      complex(-28178.111713034479,  74357.58237270877),
      complex( 14629.740252320162, -19181.808185015423),
      complex( -2870.4181610336977,  1674.1094840839117),
      complex(   132.16594124747165,   17.47674798877043)];
nullst(10)=Z10;
nullst(11)=KP10;
// M=8,N=6
Z8  =[complex(9.1616288675704660,1.5583773940642440),
      complex(8.4022529851546650,4.6911627560482549),
      complex(6.7412618178791070,7.8859496894876399),
      complex(3.6948563293866590,11.2696993000969280)];
KP8 =[complex(1003.4311174825393000,-9469.9946196260789000),
      complex(-1565.3696169904713000,4507.8810923293640000),
      complex(625.6359520664356000,-857.1892267478720000),
      complex(-63.6974526131752230,29.9056373247532870)];
nullst(8)=Z8;
nullst(9)=KP8;
// M=6,N=4
Z6  =[complex(6.4825344919041280,1.4993106204860403),
      complex(5.4692594645750200,4.5199308032962099),
      complex(3.0482060435208178,7.6517910197573519)];
KP6 =[complex(55.6824200784903400,-630.1306193716413900),
      complex(-74.3139522943279620,249.2676609831253600),
      complex(18.6315322158254460,-26.2157418820947750)];
nullst(6)=Z6;
nullst(7)=KP6;
// M=4,N=2
Z4  =[complex(3.7790199670101930,1.3801765242728430),
      complex(2.2209800329898070,4.1603914455069319)];
KP4 =[complex(2.2569587444181398,-39.6330870005017180),
      complex(-2.2569587444181418,11.1088316378758970)];
nullst(4)=Z4;
nullst(5)=KP4;
// M=2,N=0
Z2  =[complex(1.0000000000000000,1.0000000000000000)];
KP2 =[complex(0.0000000000000000,-2.0000000000000000)];
nullst(2)=Z2;
nullst(3)=KP2;
// End Zeros/Residues Library

// Get zeros and residues according to power of M and
// not very much higher.
// If more accuracy is needed "bump" the number of Zeros/Residues with
// the variable zeroPrec. Check results !
Ndash=N+zeroPrec;
if modulo(Ndash,2)==0 then
    Z  = nullst(Ndash);
    KP = nullst(Ndash+1);
else
    Z  = nullst(Ndash+1);
    KP = nullst(Ndash+2);
end
//++++++++++++++++++++++++++++++++++++++++++++
X      = zeros(15,1); // Result Vector
xa     = zeros(15,1); // Init condition answer
xnull  = zeros(15,1); // Init condition@t=0
qi     = zeros(15,length(Z));

xt  = zeros(TEND/KK+5);
vt  = zeros(TEND/KK+5);

// Compute the M' complex vectors q1,q2... in qi
for k=1:length(Z)
    s = Z(k)/h;
    qi(1,k) = A(1);
    for j=2:N
        qi(j,k) = A(j)+qi(j-1,k)/s;
    end
end
// Step 4
for t=0.0:h:TEND
    // Step 5
    xa = xa*0.0; // reset xa
    for j=1:1:length(Z) // counter for zeros\residues
        s  = Z(j)/h;
        // User supplied input function W(s)
        Ws = 1.0/s;
        // Step 6 ; Eq 22 -> sum(qi*xnull) from 1 to N-1
        summe = sum(qi(1:N-1,1)'*xnull(1:N-1));
        // calculate X(N); X=temp vector for calculus
        X(N)  = (Ws + xnull(N)-summe/s)/(s+qi(N,j));
        // End Eq 22
        // Eq 23, calculate X(N-1)..X(1)
        for k=N-1:-1:1
            X(k)  = (xnull(k) + X(k+1))/s;
        end
        // End Eq 23
        // Step 7 calculate result vector x
        xa = xa - (1/h)*real(KP(j)*X);
    end
    xnull = xa;
    //output every outKK time step for plot
    if modulo(t,KK)<0.01 then
        yt = B'*xnull; // Eq 16: v(t)=sum(Bi*xi)i=1...M
        vt = [vt,yt];  //  v(t) container
        xt = [xt,t+h]; //  time points container
    end
end

plot(xt,vt,'--b.'); // blue, dashed line + dots
xlabel("Numerical Inverse");
// end of program

//#####################################################
//#                                                   #
//#                  Scratchpad                       #
//#                                                   #
//#####################################################
//    yt    = 0.0;
//    for i=1:1:M
//        yt = yt + B(i)*xnull(i);
//    end

        // sq3=sqrt(15);
        // func=(10/sq3)*sinh((t+h)*sq3/10)*%e^(-(t+h)/2);
        // mprintf(Outformat,t+h,yt,func);

//plot(Xaxis,Yaxis,'--k',x,[(4.6*sin(0.87*x)-2*cos(0.87*x))./exp(x/2) + 2],'r+');
//plot(Xaxis,Yaxis,'--k',Xaxis,[(1-cos(Xaxis))],'r+');
//xlabel("--=Numerical Inverse, red=Inverse Laplace Function");

//        summe = 0.0;
//        for i=1:N-1
//            summe = summe + qi(i,j)*xnull(i);
//        end

// Format for console print
// Outformat = "t: %3.2f   y(t): %-08.5f  func:%-08.5f\n";
