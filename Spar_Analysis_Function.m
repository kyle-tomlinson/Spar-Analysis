function [name, PID, dataOut1, dataOut2] = Spar_Analysis_Function(dataIn)
% + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
% +
% +  SE-160A:  Aerospace Structural Analysis I
% +
% +  Project: (1) Spar Analysis
% +
% +  Title:   Spar_Analysis_Function
% +  Author:  Kyle Tomlinson
% +  PID:     A16902726
% +  Revised: 02/19/2024
% +
% +  This function is the primary analysis function for the spar analysis 
% +  program.  All of the input data is brought into the function using 
% +  "dataIn". Next all the calculations are performed.  Finally, the 
% +  calculated results are written to "dataOut1" and "dataOut2", where 
% +  these two data sets are sent to the main program (p-code) where it is 
% +  written and plotted in an Excel output file.
% +
% +  SUMMARY OF SPAR ANALYSIS
% +
% +  A) SECTION PROPERTIES
% +     A.1) Section Properties (EA EIyy, EIzz, EIyz)
% +     A.2) Torsion Constant (GJ)
% +  B) LOADS
% +     B.1) Applied Concentrated Forces, Torque and Moment
% +     B.2) Distributed Aerodynamic Loads (lift, drag, moment)
% +     B.3) Spar Root Reactions - Axial, Shear, Torque, and Moment 
% +     B.4) Distributed Spar Internal Axial, Shear, Bend Moments, Torque
% +  C) INTERNAL STRESSES (Four Cross-Section Locations)
% +     C.1) Root Spar Axial Stress (sxx) and Shear Stress (txs)
% +     C.2) Allowable Stress, and Root Margin of Safety
% +     C.3) Distributed Spar Axial Stress (sxx)
% +     C.4) Distributed Spar Shear Stress (tau)
% +  D) SPAR TIP DISPLACEMENTS, TWIST, AND BENDING SLOPES
% +     D.1) Distributed Spar X-direction (Axial) Displacement 
% +     D.2) Distributed Spar Y-direction (Drag ) Displacement 
% +     D.3) Distributed Spar Z-direction (Lift ) Displacement 
% +     D.4) Distributed Spar Twist Rotation
% +     D.5) Distributed Spar Bending Slopes (dv/dx, dw/dx)
% +
% +  Input Data
% +     dataIn:           Packed input data (38)
% +       dataIn(01):     Number of Output Plot Data Points
% +       dataIn(02):     Spar Length (inch)
% +       dataIn(03):     Mean Cross-Section Radius (inch)
% +       dataIn(04):     Mean Cross-Section Thickness (inch)
% +       dataIn(05):     Material Density (lbf/in^3)
% +       dataIn(06):     Material Young's Modulus (Msi)
% +       dataIn(07):     Material Shear Modulus (Msi)
% +       dataIn(08):     Material Yield Strength - Tension (Ksi)
% +       dataIn(09):     Material Yield Strength - Compression (Ksi)
% +       dataIn(10):     Material Yield Strength - Shear (Ksi)
% +       dataIn(11):     Material Ultimate Strength - Tension (Ksi)
% +       dataIn(12):     Material Ultimate Strength - Compression (Ksi)
% +       dataIn(13):     Material Ultimate Strength - Shear (Ksi)
% +       dataIn(14):     Safety Factor - Yield
% +       dataIn(15):     Safety Factor - Ultimate
% +       dataIn(16):     First Load Location (x/L)
% +       dataIn(17):     Concentrated Force - X Direction (lb)
% +       dataIn(18):     Concentrated Force - Y Direction (lb)
% +       dataIn(19):     Concentrated Force - Z Direction (lb) 
% +       dataIn(20):     Concentrated Torque - About X Direction (lb-in)
% +       dataIn(21):     Concentrated Moment - About Y Direction (lb-in)
% +       dataIn(22):     Concentrated Moment - About Z Direction (lb-in)
% +       dataIn(23):     Second Load Location (x/L)
% +       dataIn(24):     Concentrated Force - X Direction (lb)
% +       dataIn(25):     Concentrated Force - Y Direction (lb)
% +       dataIn(26):     Concentrated Force - Z Direction (lb) 
% +       dataIn(27):     Concentrated Torque - About X Direction (lb-in)
% +       dataIn(28):     Concentrated Moment - About Y Direction (lb-in)
% +       dataIn(29):     Concentrated Moment - About Z Direction (lb-in)
% +       dataIn(30):     Aircraft Load Factor
% +       dataIn(31):     Drag Distribution - Constant (lb/in)
% +       dataIn(32):     Drag Distribution - rth order (lb/in)
% +       dataIn(33):     Drag Distribution - polynomial order
% +       dataIn(34):     Lift Distribution - Constant (lb/in)
% +       dataIn(35):     Lift Distribution - 2nd Order (lb/in)
% +       dataIn(36):     Lift Distribution - 4th Order (lb/in)
% +       dataIn(37):     Twist Moment Distribution - Constant (lb-in/in)
% +       dataIn(38):     Twist Moment Distribution - 1st Order (lb-in/in)
% +
% +  Output Data
% +     Name:             Name of author of this analysis function           
% +     PID:              UCSD Student ID number of author
% +     dataOut1:         Packed calculated output variable data
% +       dataOut1(01):   Axial   Stiffness EA   (lb)
% +       dataOut1(02):   Bending Stiffness EIyy (lb-in^2)   
% +       dataOut1(03):   Bending Stiffness EIzz (lb-in^2)   
% +       dataOut1(04):   Bending Stiffness EIyz (lb-in^2)   
% +       dataOut1(05):   Torsion Stiffness GJ   (lb-in^2)
% +       dataOut1(06):   Root Internal Force - X Direction (lb)
% +       dataOut1(07):   Root Internal Force - Y Direction (lb)
% +       dataOut1(08):   Root Internal Force - Z Direction (lb)
% +       dataOut1(09):   Root Internal Moment - about X Direction (lb-in)
% +       dataOut1(10):   Root Internal Moment - about Y Direction (lb-in)
% +       dataOut1(11):   Root Internal Moment - about Z Direction (lb-in)
% +       dataOut1(12):   Allowable Stress - Tension (lb/in^2)
% +       dataOut1(13):   Allowable Stress - Compression (lb/in^2)
% +       dataOut1(14):   Allowable Stress - Shear (lb/in^2)
% +       dataOut1(15):   Root Axial Stress - point A (lb/in^2)
% +       dataOut1(16):   Root Axial Stress - point B (lb/in^2)
% +       dataOut1(17):   Root Axial Stress - point C (lb/in^2)
% +       dataOut1(18):   Root Axial Stress - point D (lb/in^2)
% +       dataOut1(19):   Root Shear Stress xy - point A (lb/in^2)
% +       dataOut1(20):   Root Shear Stress xy - point B (lb/in^2)
% +       dataOut1(21):   Root Shear Stress xy - point C (lb/in^2)
% +       dataOut1(22):   Root Shear Stress xy - point D (lb/in^2)
% +       dataOut1(23):   Root Shear Stress xz - point A (lb/in^2)
% +       dataOut1(24):   Root Shear Stress xz - point B (lb/in^2)
% +       dataOut1(25):   Root Shear Stress xz - point C (lb/in^2)
% +       dataOut1(26):   Root Shear Stress xz - point D (lb/in^2)
% +       dataOut1(27):   Margin of Safety - point A
% +       dataOut1(28):   Margin of Safety - point B
% +       dataOut1(29):   Margin of Safety - point C
% +       dataOut1(30):   Margin of Safety - point D
% +       dataOut1(31):   Tip Diplacement  - X Direction (inch)
% +       dataOut1(32):   Tip Diplacement  - Y Direction (inch)
% +       dataOut1(33):   Tip Diplacement  - Z Direction (inch)
% +       dataOut1(34):   Tip Twist (degree)
% +       dataOut1(35):   Tip Bending Slope (dv/dx) (inch/inch)
% +       dataOut1(36):   Tip Bending Slope (dw/dx) (inch/inch)
% + 
% +     dataOut2:         Packed calculated output plot data
% +       column( 1):     X direction coordinate (inch)
% +       column( 2):     Applied distributed drag force (lb/in)   
% +       column( 3):     Aapplied distributed lift force (lb/in)   
% +       column( 4):     Applied distributed torque (lb-in/in)   
% +       column( 5):     Internal axial force  - Vx (lb)   
% +       column( 6):     Internal shear force  - Vy (lb)   
% +       column( 7):     Internal shear force  - Vz (lb)   
% +       column( 8):     Internal axial torque - Mx (lb-in)   
% +       column( 9):     Internal bending moment - My (lb-in)
% +       column(10):     Internal bending moment - Mz (lb-in) 
% +       column(11):     Axial Stress - point A (lb/in^2) 
% +       column(12):     Axial Stress - point B (lb/in^2) 
% +       column(13):     Axial Stress - point C (lb/in^2) 
% +       column(14):     Axial Stress - point D (lb/in^2) 
% +       column(15):     Shear Stress - point A (lb/in^2) 
% +       column(16):     Shear Stress - point B (lb/in^2) 
% +       column(17):     Shear Stress - point C (lb/in^2) 
% +       column(18):     Shear Stress - point D (lb/in^2)
% +       column(19):     Displacement - X Direction (inch) 
% +       column(20):     Displacement - Y Direction (inch) 
% +       column(21):     Displacement - z Direction (inch) 
% +       column(22):     Twist (degree)
% +       column(23):     Bending Slope (dv/dx) (inch/inch)
% +       column(24):     Bending Slope (dw/dx) (inch/inch)
% +
% + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#1): Unpack Input Data Array and Wwrite User Name and PID
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (1) Unpack Input Data and Write User Name and PID')

   nplot = dataIn( 1);   % number of output plot data points
   Lo    = dataIn( 2);   % Spar Length (inch)
   Ro    = dataIn( 3);   % Mean Cross-Section Radius (inch)
   to    = 0.0285;              %dataIn( 4);   Mean Cross-Section Thickness (inch)
   rho   = dataIn( 5);   % Material Density (lbf/in^3)
   Eo    = dataIn( 6);   % Material Young's Modulus (Msi)
   Go    = dataIn( 7);   % Material Shear Modulus (Msi)
   Syt   = dataIn( 8);   % Material Yield Strength - Tension (Ksi)
   Syc   = dataIn( 9);   % Material Yield Strength - Compression (Ksi)
   Sys   = dataIn(10);   % Material Yield Strength - Shear (Ksi)
   Sut   = dataIn(11);   % Material Ultimate Strength - Tension (Ksi)
   Suc   = dataIn(12);   % Material Ultimate Strength - Compression (Ksi)
   Sus   = dataIn(13);   % Material Ultimate Strength - Shear (Ksi)
   SFy   = dataIn(14);   % Safety Factor - Yield
   SFu   = dataIn(15);   % Safety Factor - Ultimate
   x1L   = dataIn(16);   % First Load Location (x/L)
   Fx1   = dataIn(17);   % Concentrated Force - X Direction (lb)
   Fy1   = dataIn(18);   % Concentrated Force - Y Direction (lb)
   Fz1   = dataIn(19);   % Concentrated Force - Z Direction (lb)   
   Mx1   = dataIn(20);   % Concentrated Torque - About X Direction (lb-in)
   My1   = dataIn(21);   % Concentrated Moment - About Y Direction (lb-in)
   Mz1   = dataIn(22);   % Concentrated Moment - About Z Direction (lb-in)
   x2L   = dataIn(23);   % Second Load Location (x/L)
   Fx2   = dataIn(24);   % Concentrated Force - X Direction (lb)
   Fy2   = dataIn(25);   % Concentrated Force - Y Direction (lb)
   Fz2   = dataIn(26);   % Concentrated Force - Z Direction (lb)   
   Mx2   = dataIn(27);   % Concentrated Torque - About X Direction (lb-in)
   My2   = dataIn(28);   % Concentrated Moment - About Y Direction (lb-in)
   Mz2   = dataIn(29);   % Concentrated Moment - About Z Direction (lb-in)
   LF    = dataIn(30);   % Aircraft Load Factor
   py0   = dataIn(31);   % Drag Distribution - Constant (lb/in)
   pyr   = dataIn(32);   % Drag Distribution - rth order (lb/in)
   rth   = dataIn(33);   % Drag Distribution - polynomial order
   pz0   = dataIn(34);   % Lift Distribution - Constant (lb/in)
   pz2   = dataIn(35);   % Lift Distribution - 2nd Order (lb/in)
   pz4   = dataIn(36);   % Lift Distribution - 4th Order (lb/in)
   mx0   = dataIn(37);   % Twist Moment Distribution - Constant (lb-in/in)
   mx1   = dataIn(38);   % Twist Moment Distribution - 1st Order (lb-in/in)

% Define author name and PID (Write in your name and PID)    
   name  = {'Kyle Tomlinson'};
   PID   = {'A16902726'};

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#2): Calculate the Section Properties
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (2) Calculate the Section Properties')

   EA = Eo*1e6*pi*2*Ro*to;
   EIyy = Eo*1e6*(pi/4)*((Ro + to/2)^4 - (Ro - to/2)^4);
   EIzz = Eo*1e6*(pi/4)*((Ro + to/2)^4 - (Ro - to/2)^4);
   EIyz = Eo*1e6*0;
   J = (pi/2)*((Ro + to/2)^4 - (Ro - to/2)^4);
   GJ = Go*1e6*J;
  
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#3): Calculate Root Internal Stress Resultants for Applied
% .                Concentrated Forces and Applied Aerodynamic Loads
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (3) Calculate Root Stress Resultants for Applied Concentrated Loads and Aero Loads')

   Vxo = Fx1 + Fx2;

   A = pi*2*Ro*to;
   V = A*Lo;
   W = V*rho*LF;
   Pz = @(x) pz4*(x/Lo).^4 + pz2*(x/Lo).^2 + pz0;
   Lift = integral(Pz, 0, Lo);
   Vzo = Lift - W + Fz1 + Fz2;

   Py = @(x) pyr*(x/Lo).^rth + py0;
   Vyo = integral(Py, 0, Lo) + Fy1 + Fy2;

   Mxf = @(x) mx1*(x/Lo) + mx0;
   Mxo = integral(Mxf, 0, Lo) + Mx1 + Mx2;

   Mzf = @(x) py0*(x^2/2) + pyr*x^(rth+2)/((rth+2)*Lo^rth);
   Mzo = Mzf(Lo) - Mzf(0) + (Fy1*(x1L)*Lo) + (Fy2*(x2L)*Lo) + Mz1 + Mz2;

   Myf = @(x) -(((pz4*x^(6))/((6)*Lo^4)) + ((pz2*x^(4))/((4)*Lo^2)) + ((pz0*x^2)/2));
   Myo = Myf(Lo) - Myf(0) + (W*Lo)/2 - (Fz1*(x1L)*Lo) - (Fz2*(x2L)*Lo) + My1 + My2;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#4): Calculate Allowable Properties, Root Stresses and Margin
% .                of Safety
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (4) Calculate Allowable Properties, Root Stresses, and Margins of Safety')
   
    S_allow_T = min(Syt/SFy, Sut/SFu);
    S_allow_C = -1*min(abs(Syc/SFy), abs(Suc/SFu));
    S_allow_S = min(Sys/SFy, Sus/SFu);

    Sxxo_A = ((Vxo/A) + (Myo*(Ro+to/2))/(EIyy/(Eo*1e6)))*1e-3;
    Sxxo_B = ((Vxo/A) + (Mzo*(Ro+to/2))/(EIzz/(Eo*1e6)))*1e-3;
    Sxxo_C = ((Vxo/A) + (Myo*(-Ro-to/2))/(EIyy/(Eo*1e6)))*1e-3;
    Sxxo_D = ((Vxo/A) + (Mzo*(-Ro-to/2))/(EIzz/(Eo*1e6)))*1e-3;

    Txyo_A = (Vyo/(pi*Ro*to) - (Mxo * (Ro+to/2)/J))*1e-3;
    Txyo_B = 0;
    Txyo_C = (Vyo/(pi*Ro*to) + (Mxo * (Ro+to/2)/J))*1e-3;
    Txyo_D = 0;

    Txzo_A = 0;
    Txzo_B = (Vzo/(pi*Ro*to) - (Mxo * (Ro+to/2)/J))*1e-3;
    Txzo_C = 0;
    Txzo_D = (Vzo/(pi*Ro*to) + (Mxo * (Ro+to/2)/J))*1e-3;

    MS_A = (S_allow_T/(Sxxo_A^2 + 3*Txyo_A^2)^0.5) - 1;
    MS_B = (S_allow_T/(Sxxo_B^2 + 3*Txzo_B^2)^0.5) - 1;
    MS_C = (S_allow_T/(Sxxo_C^2 + 3*Txyo_C^2)^0.5) - 1;
    MS_D = (S_allow_T/(Sxxo_D^2 + 3*Txzo_D^2)^0.5) - 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#5): Calculate the Data Arrays for Plotting
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (5) Calculate the Data Arrays for Future Plotting')

   x = linspace(0, Lo, nplot);

   Vxf = @(x) Vxo;
   Vzf = @(x) Vzo -(pz0*x + (pz2*x.^3)/(3*Lo^2) + (pz4*x.^5)/(5*Lo^4) - (W*x)/Lo);
   Vyf = @(x) Vyo -((pyr*x.^7)/(7*Lo^6) + (py0*x));

   Mxfct = @(x) Mxo - (mx0*x + (mx1*x)/3);
   Mzfct = @(x) Mzo + (-Vyo*x + pyr*x.^(rth+2)/((rth+1)*(rth+2)*Lo^rth) + py0*x.^2/2);
   Myfct = @(x) Myo + Vzo*x - ((pz4/(30*Lo.^4)*x.^6 + (pz2/(12*Lo.^2))*x.^4 + (pz0/2)*x.^2)) + (W*x^2)/(2*Lo);

   for n = 1:nplot
        py(n) = Py(x(n));
        pz(n) = Pz(x(n));
        mx(n) = Mxf(x(n));
   end 

   for n = 1:nplot
       if x(n) < x1L*Lo
            Vx(n) = Vxf(x(n));
            Vy(n) = Vyf(x(n));
            Vz(n) = Vzf(x(n));
            Mx(n) = Mxfct(x(n));
            My(n) = Myfct(x(n));
            Mz(n) = Mzfct(x(n));
       end 
       if x(n) > x1L*Lo && x(n) < x2L*Lo
            Vx(n) = Vxf(x(n)) - Fx1;
            Vy(n) = Vyf(x(n)) - Fy1;
            Vz(n) = Vzf(x(n)) - Fz1;
            Mx(n) = Mxfct(x(n)) - Mx1;
            My(n) = Myfct(x(n)) + (Fz1*((x1L*Lo)+x(n))) + (Fz2*((x2L*Lo)-(x1L*Lo)-x(n)))- My1;
            Mz(n) = Mzfct(x(n)) + (Fy1*(x1L)*Lo) - Mz1;
       end 
       if x(n) > x2L*Lo
            Vx(n) = Vxf(x(n)) - Fx1 - Fx2;
            Vy(n) = Vyf(x(n)) - Fy1 - Fy2;
            Vz(n) = Vzf(x(n)) - Fz1 - Fz2;
            Mx(n) = Mxfct(x(n)) - Mx1 - Mx2;
            My(n) = Myfct(x(n)) + (Fz2*(x2L)*Lo) - My1;
            Mz(n) = Mzfct(x(n)) + (Fz1*((x2L*Lo)+x(n))) + (Fz2*x(n))- My1 - Mz2; 
       else
            Vx(n) = Vxf(x(n));
            Vy(n) = Vyf(x(n));
            Vz(n) = Vzf(x(n));
            Mx(n) = Mxfct(x(n));
            My(n) = Myfct(x(n));
            Mz(n) = Mzfct(x(n));
       end 
   end 
   ((Vxo/A) + (Myo*(Ro+to/2))/(EIyy/(Eo*1e6)))*1e-3;
   SxxAf = @(x) My(x)*(Ro+to/2)/(EIyy/(Eo*1e6))*1e-3;
   SxxBf = @(x) Mz(x)*(Ro+to/2)/(EIzz/(Eo*1e6))*1e-3;
   SxxCf = @(x) My(x)*(-Ro-to/2)/(EIyy/(Eo*1e6))*1e-3;
   SxxDf = @(x) Mz(x)*(-Ro-to/2)/(EIzz/(Eo*1e6))*1e-3;

   TauAf = @(x) ((Vy(x)/(pi*Ro*to)) - ((Mx(x) * (Ro+to/2)/J)))*1e-3;
   TauBf = @(x) ((Vz(x)/(pi*Ro*to)) - ((Mx(x) * (Ro+to/2)/J)))*1e-3;
   TauCf = @(x) ((Vy(x)/(pi*Ro*to)) + ((Mx(x) * (Ro+to/2)/J)))*1e-3;
   TauDf = @(x) ((Vz(x)/(pi*Ro*to)) + ((Mx(x) * (Ro+to/2)/J)))*1e-3;

   for n = 1:nplot
        Sxx_A(n) = SxxAf(n);
        Sxx_B(n) = SxxBf(n);
        Sxx_C(n) = SxxCf(n);
        Sxx_D(n) = SxxDf(n);
        Tau_A(n) = TauAf(n);
        Tau_B(n) = TauBf(n);
        Tau_C(n) = TauCf(n);
        Tau_D(n) = TauDf(n);
   end 

   dvdxf = @(x) (Mzo*x + (-Vyo*x^2)/2 + (pyr*x.^(rth+3))/((rth+2)*(rth+1)*(rth+3)*Lo^rth) + py0*x.^3/6)/(EIyy);
   dwdxf = @(x) (Myo*x + (Vzo*x^2)/2 - (pz4*x^7)/(210*Lo^4) - (pz2*x^5)/(60*Lo^2) - (pz0*x^3)/6 + (W*x^3)/(Lo*6))/(-EIyy);
   dispxf = @(x) ((Fx1 + Fx2)*Lo)/EA;
   dispyf = @(x) ((Mzo*x^2)/2 + (-Vyo*x^3)/6 + (pyr*x.^(rth+4))/((rth+4)*(rth+2)*(rth+1)*(rth+3)*Lo^rth) + py0*x.^4/24)/(EIyy);
   dispzf = @(x) ((Myo*x^2)/2 + (Vzo*x^3)/6 - (pz4*x^8)/(1680*Lo^4) - (pz2*x^6)/(360*Lo^2) - (pz0*x^4)/24 + (W*x^4)/(Lo*24))/(-EIyy);
   twistf = @(x) ((Mxo*x) - (mx1*x^2)/6 - (mx0*x^2)/2)*(1/GJ) * (180/pi);
   
   for n = 1:nplot
        Disp_X(n) = dispxf(x(n));
        Disp_Y(n) = dispyf(x(n));
        Disp_Z(n) = dispzf(x(n));
        DvDx(n) = dvdxf(x(n));
        DwDx(n) = dwdxf(x(n));
        Twist(n) = twistf(x(n));
   end 

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#6): Pack Calculated Data into the "dataOut1" Array size: (36)
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (6) Pack the Calculated Data into Array: dataOut1')

     dataOut1(01) = EA;              % Axial   Stiffness (lb)
     dataOut1(02) = EIyy;            % Bending Stiffness (lb-in^2)   
     dataOut1(03) = EIzz;            % Bending Stiffness (lb-in^2)   
     dataOut1(04) = EIyz;            % Bending Stiffness (lb-in^2)   
     dataOut1(05) = GJ;              % Torsion Stiffness (lb-in^2)
     dataOut1(06) = Vxo;             % Root Internal Force - X Direction (lb)
     dataOut1(07) = Vyo;             % Root Internal Force - Y Direction (lb)
     dataOut1(08) = Vzo;             % Root Internal Force - Z Direction (lb)
     dataOut1(09) = Mxo;             % Root Internal Moment - about X Direction (lb-in)
     dataOut1(10) = Myo;             % Root Internal Moment - about Y Direction (lb-in)
     dataOut1(11) = Mzo;             % Root Internal Moment - about Z Direction (lb-in)
     dataOut1(12) = S_allow_T;       % Allowable Stress - Tension (Ksi)
     dataOut1(13) = S_allow_C;       % Allowable Stress - Compression (Ksi)
     dataOut1(14) = S_allow_S;       % Allowable Stress - Shear (Ksi)
     dataOut1(15) = Sxxo_A;          % Root Axial Stress - point A (Ksi)
     dataOut1(16) = Sxxo_B;          % Root Axial Stress - point B (Ksi)
     dataOut1(17) = Sxxo_C;          % Root Axial Stress - point C (Ksi)
     dataOut1(18) = Sxxo_D;          % Root Axial Stress - point D (Ksi)
     dataOut1(19) = Txyo_A;          % Root Shear Stress xy - point A (ksi)   
     dataOut1(20) = Txyo_B;          % Root Shear Stress xy - point B (ksi)
     dataOut1(21) = Txyo_C;          % Root Shear Stress xy - point C (ksi)
     dataOut1(22) = Txyo_D;          % Root Shear Stress xy - point D (ksi)
     dataOut1(23) = Txzo_A;          % Root Shear Stress xz - point A (ksi)
     dataOut1(24) = Txzo_B;          % Root Shear Stress xz - point B (ksi)
     dataOut1(25) = Txzo_C;          % Root Shear Stress xz - point C (ksi)
     dataOut1(26) = Txzo_D;          % Root Shear Stress xz - point D (ksi)
     dataOut1(27) = MS_A;            % Margin of Safety - point A
     dataOut1(28) = MS_B;            % Margin of Safety - point B
     dataOut1(29) = MS_C;            % Margin of Safety - point C
     dataOut1(30) = MS_D;            % Margin of Safety - point D
     dataOut1(31) = Disp_X(nplot);   % Tip Diplacement - Y Direction (inch)
     dataOut1(32) = Disp_Y(nplot);   % Tip Diplacement - Y Direction (inch)
     dataOut1(33) = Disp_Z(nplot);   % Tip Diplacement - Z Direction (inch)
     dataOut1(34) = Twist(nplot);    % Tip Twist (degree)
     dataOut1(35) = DvDx(nplot);     % Tip Bending Slope (dv/dx) (inch/inch)
     dataOut1(36) = DwDx(nplot);     % Tip Bending Slope (dw/dx) (inch/inch)

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% .  SECTION (#7): Pack the plot data arrays into "dataOut2" 
% .                matrix size: (nplot,23)  
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   disp('     (7) Pack the Calculated Plot Data into Array: dataOut2')

   for i = 1:nplot
     dataOut2(i, 1) = x(i);          % x-location (inch)
     dataOut2(i, 2) = py(i);         % drag force (lb/in)
     dataOut2(i, 3) = pz(i);         % lift force (lb/in)
     dataOut2(i, 4) = mx(i);         % distributed torque (lb-in/in)
     dataOut2(i, 5) = Vx(i);         % Internal axial force - Vx (lb)
     dataOut2(i, 6) = Vy(i);         % Internal shear force - Vy (lb)
     dataOut2(i, 7) = Vz(i);         % Internal shear force - Vz (lb)
     dataOut2(i, 8) = Mx(i);         % Internal axial force - Mx (lb-in)
     dataOut2(i, 9) = My(i);         % Internal shear force - My (lb-in)
     dataOut2(i,10) = Mz(i);         % Internal shear force - Mz (lb-in)
     dataOut2(i,11) = Sxx_A(i);      % Axial Stress - Point A (Ksi)
     dataOut2(i,12) = Sxx_B(i);      % Axial Stress - Point B (Ksi)
     dataOut2(i,13) = Sxx_C(i);      % Axial Stress - Point C (Ksi)
     dataOut2(i,14) = Sxx_D(i);      % Axial Stress - Point D (Ksi)
     dataOut2(i,15) = Tau_A(i);      % Shear Stress - Point A (Ksi)
     dataOut2(i,16) = Tau_B(i);      % Shear Stress - Point B (Ksi)
     dataOut2(i,17) = Tau_C(i);      % Shear Stress - Point C (Ksi)
     dataOut2(i,18) = Tau_D(i);      % Shear Stress - Point D (Ksi)
     dataOut2(i,19) = Disp_X(i);     % Displacement - X Direction (inch)
     dataOut2(i,20) = Disp_Y(i);     % Displacement - Y Direction (inch)
     dataOut2(i,21) = Disp_Z(i);     % Displacement - Z Direction (inch)
     dataOut2(i,22) = Twist(i);      % Spar Twist (degree)
     dataOut2(i,23) = DvDx(i);       % Bending Slope (dv/dx) (inch/inch)
     dataOut2(i,24) = DwDx(i);       % Bending Slope (dw/dx) (inch/inch)
   
   end
   
end

%  End of Function: SE160A_1_sparAnalysis
%  ------------------------------------------------------------------------

