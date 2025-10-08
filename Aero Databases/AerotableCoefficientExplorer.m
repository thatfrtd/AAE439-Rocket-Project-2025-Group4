aerotable_name = "hitech_Aerotable_6000-40ALT_40PHI";

%% Load Aerotable
if isfile(string(aerotable_name + ".mat"))
    aerotable = load(aerotable_name).aerotable;
else
    error("Could not find requested aerotable.")
end
%%
% 
% aerotable = aerotables(3).Aerotable;

%% Extract Values from Aerotable
ALPHA = aerotable.ALPHA_axis3;
PHI = aerotable.PHI_axis2;
MACH = aerotable.MACH_axis4;
ALT = aerotable.ALT_axis1;

%% Plot
ALT_index = 1;
signed_log = false;

f = figure;
ax = axes('Parent',f,'position',[0.13 0.36  0.77 0.54]);
view(3);

l = uicontrol(f, 'Style', 'checkbox', 'String', "Signed Log Scale");
l.Position = [20 35 120 20];
l.Value = false;

c = uicontrol(f,'Style','popupmenu');
c.Position = [20 55 60 20];
c.String = {'CA',
            'CY',
            'CN',
            'CLL',
            'CM',
            'CLN',
            'CMQ',
            'CMAD',
            'CLP',
            'CLR',
            'CNP',
            'CNR',
            'CAQ',
            'CNQ',
            'CNAD',
            'CYR',
            'CYP',
            'CLB',
            'CYB',
            'CNB',
            'XCP',
            'CD',
            'CL',
            'CN'};

b = uicontrol('Parent',f,'Style','slider','Position',[81,84,419,23],...
              'value',ALT_index, 'min',1, 'max', numel(ALT));

plotALPHA_PHI_MACH(ALPHA, PHI, MACH, aerotable, ALT, b, c, l, ax)

bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,80,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,80,50,23],...
                'String',sprintf("%g", max(ALT)),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,55,100,23],...
                'String','Altitude [feet]','BackgroundColor',bgcolor);

l.Callback = @(es,ed) plotALPHA_PHI_MACH(ALPHA, PHI, MACH, aerotable, ALT, b, c, l, ax);
c.Callback = @(es,ed) plotALPHA_PHI_MACH(ALPHA, PHI, MACH, aerotable, ALT, b, c, l, ax);
b.Callback = @(es,ed) plotALPHA_PHI_MACH(ALPHA, PHI, MACH, aerotable, ALT, b, c, l, ax); 


function [] = plotALPHA_PHI_MACH(ALPHA, PHI, MACH, aerotable, ALT, b, c, l, ax)
    cla(ax)

    coefficient_name = c.String{c.Value};
    V = aerotable.(coefficient_name + "_sustainer");
    
    ALT_index = round(b.Value);

    signed_log_scale = l.Value;

    layers = 5;

    [X,Y,Z] = meshgrid(ALPHA, PHI, MACH);
    V = squeeze(V(ALT_index,:,:,:));
    if signed_log_scale
        V = log10(abs(V)+1) .* sign(V);
    end    
    
    s = isosurface(X,Y,Z,V,0);
    p = patch(s);
    isonormals(X,Y,Z,V,p)

    set(p,'FaceColor',[0.5 1 0.5]);  
    set(p,'EdgeColor','none');
    camlight;
    lighting gouraud;
    alpha 0.5;

    xslice = ALPHA;   
    yslice = [];
    zslice = [];
    contourslice(X,Y,Z,V,xslice,yslice,zslice, layers);
    
    colorbar(ax)
    grid on

    xlabel("ALPHA [degrees]")
    ylabel("PHI [degrees]")
    zlabel("MACH")
    title(sprintf("%s vs ALPHA, PHI, & MACH at %g feet ALT", coefficient_name, ALT(ALT_index)))
    subtitle(sprintf("Transparent Green Surface is the %s = 0 Isosurface", coefficient_name))
end