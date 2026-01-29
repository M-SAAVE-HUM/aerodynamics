global xpoints ypoints zpoints %Lord forgive me for using globals

%%%%%%% MAKING UI %%%%%%%
UIFigure = uifigure('name','Airfoil Point Generator','Visible','off');
set(UIFigure, 'position',[100 100 656 378]);

%%%%%%% MAKING PLOT %%%%%%%
UIAxes = uiaxes(UIFigure);
UIAxes.XTick = [];
UIAxes.YTick = [];
UIAxes.ZTick = [];
UIAxes.Position = [282 3 375 375];

%%%%%%% MAKING SPACING BUTTON %%%%%%%
SpacingTypeButtonGroup = uibuttongroup(UIFigure);
SpacingTypeButtonGroup.Title = 'Spacing Type';
SpacingTypeButtonGroup.Position = [23 67 123 73];

CosineButton = uiradiobutton(SpacingTypeButtonGroup);
CosineButton.Text = 'Cosine';
CosineButton.Position = [11 27 59 22];
CosineButton.Value = true;

LinearButton = uiradiobutton(SpacingTypeButtonGroup);
LinearButton.Text = 'Linear';
LinearButton.Position = [11 5 55 22];

%%%%%%% MAKING TRAILING EDGE BUTTON %%%%%%%
TrailingEdgeTypeButtonGroup = uibuttongroup(UIFigure);
TrailingEdgeTypeButtonGroup.Title = 'Trailing Edge Type';
TrailingEdgeTypeButtonGroup.Position = [158 67 123 73];

ClosedButton = uiradiobutton(TrailingEdgeTypeButtonGroup);
ClosedButton.Text = 'Closed';
ClosedButton.Position = [11 27 60 22];
ClosedButton.Value = true;

ExactFormButton = uiradiobutton(TrailingEdgeTypeButtonGroup);
ExactFormButton.Text = 'Exact Form';
ExactFormButton.Position = [11 5 84 22];

%%%%%%% MAKING AIRFOIL INPUT %%%%%%%
EnterNACA4DigitNACAAirfoilLabel = uilabel(UIFigure);
EnterNACA4DigitNACAAirfoilLabel.Position = [19 244 184 53];
EnterNACA4DigitNACAAirfoilLabel.Text = 'Enter NACA 4-Digit NACA Airfoil:';

EnterNACA4DigitNACAAirfoilEditField = uieditfield(UIFigure, 'numeric');
EnterNACA4DigitNACAAirfoilEditField.Limits = [1000 9999];
EnterNACA4DigitNACAAirfoilEditField.RoundFractionalValues = 'on';
EnterNACA4DigitNACAAirfoilEditField.HorizontalAlignment = 'left';
EnterNACA4DigitNACAAirfoilEditField.Placeholder = '2412';
EnterNACA4DigitNACAAirfoilEditField.Position = [210 256 65 27];
EnterNACA4DigitNACAAirfoilEditField.Value = 2412;

%%%%%%% MAKING GENERATE BUTTON %%%%%%%
GenerateButton = uibutton(UIFigure, 'push');
GenerateButton.BackgroundColor = [0.3922 0.8314 0.0745];
GenerateButton.Position = [35 24 99 26];
GenerateButton.Text = 'Generate';

%%%%%%% MAKING EXPORT BUTTON %%%%%%%
ExportPointsButton = uibutton(UIFigure, 'push');
ExportPointsButton.BackgroundColor = [0.2392 0.6039 0.851];
ExportPointsButton.Position = [170 24 99 26];
ExportPointsButton.Text = 'Export Points';

%%%%%%% MAKING NUM POINTS INPUT %%%%%%%
EnterNumberofPointsLabel = uilabel(UIFigure);
EnterNumberofPointsLabel.Position = [19 192 197 53];
EnterNumberofPointsLabel.Text = 'Enter Number of Points (even):';

EnterNumberofPointsevenEditField = uieditfield(UIFigure, 'numeric');
EnterNumberofPointsevenEditField.Limits = [0 Inf];
EnterNumberofPointsevenEditField.HorizontalAlignment = 'left';
EnterNumberofPointsevenEditField.Placeholder = '50';
EnterNumberofPointsevenEditField.Position = [210 204 65 27];
EnterNumberofPointsevenEditField.Value = 50;

%%%%%%% MAKING CHORD LENGTH INPUT %%%%%%%
ChordLengthmmEditFieldLabel = uilabel(UIFigure);
ChordLengthmmEditFieldLabel.Position = [19 140 197 53];
ChordLengthmmEditFieldLabel.Text = 'Chord Length (mm)';

ChordLengthmmEditField = uieditfield(UIFigure, 'numeric');
ChordLengthmmEditField.Limits = [0 Inf];
ChordLengthmmEditField.HorizontalAlignment = 'left';
ChordLengthmmEditField.Placeholder = '100';
ChordLengthmmEditField.Position = [210 152 65 27];
ChordLengthmmEditField.Value = 100;

%%%%%%% MAKING TITLE BUTTON %%%%%%%
AirfoilPointGeneratorLewieLawrence15May2024Label = uilabel(UIFigure);
AirfoilPointGeneratorLewieLawrence15May2024Label.HorizontalAlignment = 'center';
AirfoilPointGeneratorLewieLawrence15May2024Label.FontName = 'Comic Sans MS';
AirfoilPointGeneratorLewieLawrence15May2024Label.Position = [80 304 134 44];
AirfoilPointGeneratorLewieLawrence15May2024Label.Text = ...
    {'Airfoil Point Generator'; 'Lewie Lawrence'; '15 May 2024'};

%%%%%%% ADDING FUNCTIONS %%%%%%%
GenerateButton.ButtonPushedFcn = @(src,event) GenerateButtonPushed( ...
    EnterNACA4DigitNACAAirfoilEditField, EnterNumberofPointsevenEditField, ...
    CosineButton, ClosedButton, UIAxes);
ExportPointsButton.ButtonPushedFcn = @(src, event) ExportPointsButtonPushed( ...
    EnterNumberofPointsevenEditField, ChordLengthmmEditField);

%%%%%%% SHOW FIGURE %%%%%%%
UIFigure.Visible = 'on';
GenerateButtonPushed( ...
    EnterNACA4DigitNACAAirfoilEditField, EnterNumberofPointsevenEditField, ...
    CosineButton, ClosedButton, UIAxes)
setMousePointer(UIFigure,GenerateButton,ExportPointsButton);


function GenerateButtonPushed(EnterNACA4DigitNACAAirfoilEditField, ...
    EnterNumberofPointsevenEditField, CosineButton, ClosedButton, UIAxes, event)
global xpoints ypoints zpoints
    %%%%%%% GATHERING USER INPUT %%%%%%%
    airfoil = EnterNACA4DigitNACAAirfoilEditField.Value;
    if(airfoil < 1000 || airfoil > 9999)
        disp("Invalid Airfoil");
        return;
    end
    thickness = mod(airfoil,100);
    camberdistance = mod(airfoil-thickness,1000)/100;
    camber = (airfoil-thickness-camberdistance*100)/1000;
    
    %%%%%%% GENERATING AIRFOIL %%%%%%%
    numPoints = EnterNumberofPointsevenEditField.Value;
    x = linspace(0,1,numPoints/2+1);
    if (CosineButton.Value == true)
        x = 0.5-0.5*cos(x*pi); % Cosine Spacing
    end
    yc = zeros(1,numPoints/2+1);
    ycprime = zeros(1,numPoints/2+1);
    m = camber/100;
    p = camberdistance/10;
    t = thickness/100;
    for i = 1:1:numPoints/2+1
        if(x(i) <= p)
            yc(i) = m/p^2*(2*p*x(i)-x(i)^2);
            ycprime(i) = 2*m/p^2*(p-x(i));
        else
            yc(i) = m/(1-p)^2*((1-2*p)+2*p*x(i)-x(i)^2);
            ycprime(i) = 2*m/(1-p)^2*(p-x(i));
        end
    end
    if (ClosedButton.Value == true)
        yt = 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4); % Closed Airfoil
    else
        yt = 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4); % Closed Airfoil
    end
    theta = atan(ycprime);
    xU = x - yt.*sin(theta);
    xL = x + yt.*sin(theta);
    yU = yc + yt.*cos(theta);
    yL = yc - yt.*cos(theta);
    
    %%%%%%% PLOTTING AIRFOIL %%%%%%%
    plot(UIAxes,x,yc,"-r");
    hold(UIAxes,"on");
    plot(UIAxes,xU,yU,"-ob");
    plot(UIAxes,xL,yL,"-ob");
    xlim(UIAxes,[-0.1,1.1]);
    ylim(UIAxes,[-0.5,0.5]);
    hold(UIAxes,"off");

    %%%%%%% SAVING AIRFOIL POINTS %%%%%%%
    xpoints = [xU,xL(2:end)];
    ypoints = [yU,yL(2:end)];
    zpoints = zeros(1,numPoints+1);
    return 
end

function ExportPointsButtonPushed(EnterNumberofPointsevenEditField, ChordLengthmmEditField, event)
global xpoints ypoints zpoints
    file = fopen("AirfoilPoints.txt", "wt");
    numpoints = EnterNumberofPointsevenEditField.Value + 1;
    chord = ChordLengthmmEditField.Value;
    for i = 1:numpoints
        fprintf( file, '%f %f %f\n', xpoints(i)*chord, ypoints(i)*chord, zpoints(i)*chord);
    end
    fclose(file);
end

function setMousePointer(fig,GenerateButton,ExportPointsButton)
fig.WindowButtonMotionFcn = @mouseMoved;

btnX = GenerateButton.Position(1);
btnY = GenerateButton.Position(2);
btnWidth = GenerateButton.Position(3);
btnHeight = GenerateButton.Position(4);

btn2X = ExportPointsButton.Position(1);
btn2Y = ExportPointsButton.Position(2);
btn2Width = ExportPointsButton.Position(3);
btn2Height = ExportPointsButton.Position(4);

    function mouseMoved(src,event)
        mousePos = fig.CurrentPoint;
        if ((mousePos(1) >= btnX) && (mousePos(1) <= btnX + btnWidth) ...
                && (mousePos(2) >= btnY) && (mousePos(2) <= btnY + btnHeight) ...
                || (mousePos(1) >= btn2X) && (mousePos(1) <= btn2X + btn2Width) ...
                && (mousePos(2) >= btn2Y) && (mousePos(2) <= btn2Y + btn2Height))
              fig.Pointer = 'hand';
        else
              fig.Pointer = 'arrow';
        end
    end
end