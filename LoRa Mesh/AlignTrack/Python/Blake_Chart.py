from MatlabRadar import *

# def blakechart(vcp, vcpangles, **args):
#     '''
#     %blakechart Range-angle-height (Blake) chart
#     %   blakechart(VCP,VCPANGLES) plots the vertical coverage pattern, VCP, for
#     %   a radar system on a Blake chart, also known as a range-height-angle
#     %   chart. The vertical coverage pattern is a detection or
#     %   constant-signal-level contour. Normal atmospheric refraction is taken
#     %   into account through the use of the CRPL exponential reference
#     %   atmosphere with a refractivity of 313 N-units and a refraction exponent
#     %   (decay constant) of 0.143859 1/km. Scattering and ducting are assumed
#     %   to be negligible.
#     %   
#     %   The range in the range-height-angle chart is the propagated range and
#     %   the height is relative to the origin of the ray. It is assumed that the
#     %   antenna height is not at an appreciable height above ground level (<
#     %   1000 ft or about 305 m).
#     %
#     %   VCP can be a matrix whose columns represent individual vertical
#     %   coverage patterns. VCPANGLES is a column vector whose number of rows is
#     %   the same as number of rows of VCP. Each entry in VCPANGLES specifies
#     %   the elevation angle (in degrees) at which the vertical coverage pattern
#     %   is measured.
#     %
#     %   blakechart(VCP,VCPANGLES,RMAX,HMAX) specifies the range limit RMAX and
#     %   height limit HMAX for the Blake chart, respectively. The units of RMAX
#     %   and HMAX are determined by the values of the RangeUnit and HeightUnit
#     %   parameters. The default unit for both is 'km'.
#     % 
#     %   blakechart(...,'RangeUnit',RNGUNIT,'HeightUnit',HTUNIT) specifies the
#     %   range and height units in RNGUNIT and HTUNIT as one of 'km' | 'm' |
#     %   'mi' | 'nmi' | 'ft' | 'kft', where the default values for both are
#     %   'km'.
#     % 
#     %   blakechart(...,'ScalePower',SPOW) specifies the range and height axis
#     %   scale power as a scalar between 0 and 1. The default value for SPOW is
#     %   1/4.
#     % 
#     %   blakechart(...,'SurfaceRefractivity',NS,'RefractionExponent',REXP)
#     %   specifies the surface refractivity NS (in N-Units) as a non-negative
#     %   scalar and the refraction exponent factor (decay constant) REXP (in
#     %   1/km) as a non-negative scalar for the CRPL exponential reference
#     %   atmosphere model. The default value of NS is 313 and the default value
#     %   of REXP is 0.143859.
#     %
#     %   The atmospheric refraction model is CRPL exponential reference model
#     %   given by
#     %        n(h) = 1 + NS*1e-6*exp(-REXP*h),
#     %   where h is the height in kilometers and n is the index of refraction.
#     %
#     %   blakechart(...,'AntennaHeight',ANHT) specifies the antenna height. The
#     %   unit of height is specified by HTUNIT. When an antenna height is
#     %   provided, the height in the Blake chart is the height above ground
#     %   level. Otherwise, the height in the Blake chart is relative to the
#     %   origin of the ray, and it is assumed that the antenna height is not at
#     %   an appreciable height above ground level (< 1000 ft or about 305 m).
#     %   Default antenna height is 0.
#     %
#     %   blakechart(...,'FaceColor',FC) specifies the face color of the vertical
#     %   coverage pattern patch. FC can be specified by the color name (e.g.,
#     %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
#     %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
#     %   of colors must match the number of columns in the VCP input. FC can
#     %   also be set to 'none' if no patch fill is desired. FC defaults to the
#     %   default MATLAB color order. See MATLAB ColorSpec for additional
#     %   information.
#     %
#     %   blakechart(...,'EdgeColor',EC) specifies the edge color of the vertical
#     %   coverage pattern patch. EC can be specified by the color name (e.g.,
#     %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
#     %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
#     %   of colors must match the number of columns in the VCP input. EC can
#     %   also be set to 'none' if no edge color is desired. EC defaults to the
#     %   default MATLAB color order. See MATLAB ColorSpec for additional
#     %   information.
#     %
#     %   blakechart(...,'Parent',HAX) specifies the plot axes, HAX. The default
#     %   axes are in the current figure.
#     % 
#     %   % Example:
#     %   %   Plot the radar vertical coverage pattern assuming the antenna has a
#     %   %   sinc pattern. The frequency is 100 MHz, the antenna height is 20
#     %   %   feet, and the range is 100 nautical miles. Assume the surface is 
#     %   %   smooth, the antenna is not tilted, and the transmitted polarization
#     %   %   is horizontal.
#     %
#     %   patAng    = linspace(-90,90,361)';
#     %   pat_u     = 1.39157/sind(90/2)*sind(patAng);
#     %   pat       = sinc(pat_u/pi);
#     % 
#     %   freq      = 100e6;
#     %   anht      = 20;
#     %   rfs       = 100;
#     %   tiltAng   = 0;
#     %   [vcdNmi, vcdAng] = radarvcd(freq,rfs,anht, ...
#     %       'RangeUnit','nmi','HeightUnit','ft', ...
#     %       'AntennaPattern',pat,'PatternAngles',patAng,'TiltAngle',tiltAng);
#     %
#     %   blakechart(vcdNmi, vcdAng, ...
#     %       'RangeUnit','nmi','HeightUnit','ft')
#     %   
#     %   See also radar, radarvcd, refractionexp.
    
#     %   Copyright 2012-2023 The MathWorks, Inc.
    
#     %   References:
#     %
#     %   [1] Blake, L.V. "Radio Ray (Radar) Range-Height-Angle Charts." Naval
#     %       Research Laboratory, NRL Report 6650, Jan. 22, 1968.
#     %   [2] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
#     %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
#     %       Series), No. 1, Jan. 1968, pp. 85-92.
#     %   [3] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
#     %       Atmosphere. Washington, DC: US Gov. Print. Off., 1959.
#     '''
    
#     # 1. Convert to Polar Coordinates
#     # Plotting elevation (theta) vs range (r)
#     r = vcp
#     theta = np.radians(vcpangles)

#     # Bookend data with zeros if not already
#     # if not (r[0] == 0) or not (r[-1] == 0) or not (theta[0] == 0) or not (theta[1] == 0):
#     #     r = np.insert(r, 0, 0)
#     #     r = np.append(r, 0)
#     #     theta = np.insert(theta, 0, 0)
#     #     theta = np.append(theta, 0)

#     r[0] = 0
#     r[1] = 0
    
#     # 2. Create Plot
#     fig, ax = plt.subplots(figsize = (8, 5), subplot_kw={'projection': 'polar'})
#     ax.plot(theta, r)
    
#     # 3. Configure Chart
#     ax.set_theta_zero_location('E')  # Set 0 degrees to North (top)
#     ax.set_theta_direction(1)       # Clockwise
#     ax.set_rmax(max(r)*1.1)
#     ax.set_rmin(0)
#     ax.set_thetamin(0)
#     ax.set_thetamax(90)
#     # ax.set_rscale('log')
#     ax.set_title("Blake Chart Representation", va='bottom')
#     ax.fill(theta, r, color='C0', alpha=0.5)
    
#     plt.show()

#     return args

def blakechart(
    vcp, # Vertical coverage pattern
    vcpangles, # Vertical coverage pattern angles
    rmax = None, # Maximum range of plot
    hmax = None, # Maximum height of plot
    rngunit = 'km', # RangeUnit
    htunit = 'km', # HeightUnit
    spow = 0.25, # ScalePower
    Ns = 313, # SurfaceRefractivity
    expo = 0.143859, # RefractionExponent
    anht = 0, # AntennaHeight
    fc = 'blue', # FaceColor
    ec = 'blue', # EdgeColor
    hAxes = None, # Parent axes
):
    '''
    %blakechart Range-angle-height (Blake) chart
    %   blakechart(VCP,VCPANGLES) plots the vertical coverage pattern, VCP, for
    %   a radar system on a Blake chart, also known as a range-height-angle
    %   chart. The vertical coverage pattern is a detection or
    %   constant-signal-level contour. Normal atmospheric refraction is taken
    %   into account through the use of the CRPL exponential reference
    %   atmosphere with a refractivity of 313 N-units and a refraction exponent
    %   (decay constant) of 0.143859 1/km. Scattering and ducting are assumed
    %   to be negligible.
    %   
    %   The range in the range-height-angle chart is the propagated range and
    %   the height is relative to the origin of the ray. It is assumed that the
    %   antenna height is not at an appreciable height above ground level (<
    %   1000 ft or about 305 m).
    %
    %   VCP can be a matrix whose columns represent individual vertical
    %   coverage patterns. VCPANGLES is a column vector whose number of rows is
    %   the same as number of rows of VCP. Each entry in VCPANGLES specifies
    %   the elevation angle (in degrees) at which the vertical coverage pattern
    %   is measured.
    %
    %   blakechart(VCP,VCPANGLES,RMAX,HMAX) specifies the range limit RMAX and
    %   height limit HMAX for the Blake chart, respectively. The units of RMAX
    %   and HMAX are determined by the values of the RangeUnit and HeightUnit
    %   parameters. The default unit for both is 'km'.
    % 
    %   blakechart(...,'RangeUnit',RNGUNIT,'HeightUnit',HTUNIT) specifies the
    %   range and height units in RNGUNIT and HTUNIT as one of 'km' | 'm' |
    %   'mi' | 'nmi' | 'ft' | 'kft', where the default values for both are
    %   'km'.
    % 
    %   blakechart(...,'ScalePower',SPOW) specifies the range and height axis
    %   scale power as a scalar between 0 and 1. The default value for SPOW is
    %   1/4.
    % 
    %   blakechart(...,'SurfaceRefractivity',NS,'RefractionExponent',REXP)
    %   specifies the surface refractivity NS (in N-Units) as a non-negative
    %   scalar and the refraction exponent factor (decay constant) REXP (in
    %   1/km) as a non-negative scalar for the CRPL exponential reference
    %   atmosphere model. The default value of NS is 313 and the default value
    %   of REXP is 0.143859.
    %
    %   The atmospheric refraction model is CRPL exponential reference model
    %   given by
    %        n(h) = 1 + NS*1e-6*exp(-REXP*h),
    %   where h is the height in kilometers and n is the index of refraction.
    %
    %   blakechart(...,'AntennaHeight',ANHT) specifies the antenna height. The
    %   unit of height is specified by HTUNIT. When an antenna height is
    %   provided, the height in the Blake chart is the height above ground
    %   level. Otherwise, the height in the Blake chart is relative to the
    %   origin of the ray, and it is assumed that the antenna height is not at
    %   an appreciable height above ground level (< 1000 ft or about 305 m).
    %   Default antenna height is 0.
    %
    %   blakechart(...,'FaceColor',FC) specifies the face color of the vertical
    %   coverage pattern patch. FC can be specified by the color name (e.g.,
    %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
    %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
    %   of colors must match the number of columns in the VCP input. FC can
    %   also be set to 'none' if no patch fill is desired. FC defaults to the
    %   default MATLAB color order. See MATLAB ColorSpec for additional
    %   information.
    %
    %   blakechart(...,'EdgeColor',EC) specifies the edge color of the vertical
    %   coverage pattern patch. EC can be specified by the color name (e.g.,
    %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
    %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
    %   of colors must match the number of columns in the VCP input. EC can
    %   also be set to 'none' if no edge color is desired. EC defaults to the
    %   default MATLAB color order. See MATLAB ColorSpec for additional
    %   information.
    %
    %   blakechart(...,'Parent',HAX) specifies the plot axes, HAX. The default
    %   axes are in the current figure.
    % 
    %   % Example:
    %   %   Plot the radar vertical coverage pattern assuming the antenna has a
    %   %   sinc pattern. The frequency is 100 MHz, the antenna height is 20
    %   %   feet, and the range is 100 nautical miles. Assume the surface is 
    %   %   smooth, the antenna is not tilted, and the transmitted polarization
    %   %   is horizontal.
    %
    %   patAng    = linspace(-90,90,361)';
    %   pat_u     = 1.39157/sind(90/2)*sind(patAng);
    %   pat       = sinc(pat_u/pi);
    % 
    %   freq      = 100e6;
    %   anht      = 20;
    %   rfs       = 100;
    %   tiltAng   = 0;
    %   [vcdNmi, vcdAng] = radarvcd(freq,rfs,anht, ...
    %       'RangeUnit','nmi','HeightUnit','ft', ...
    %       'AntennaPattern',pat,'PatternAngles',patAng,'TiltAngle',tiltAng);
    %
    %   blakechart(vcdNmi, vcdAng, ...
    %       'RangeUnit','nmi','HeightUnit','ft')
    %   
    %   See also radar, radarvcd, refractionexp.
    
    %   Copyright 2012-2023 The MathWorks, Inc.
    
    %   References:
    %
    %   [1] Blake, L.V. "Radio Ray (Radar) Range-Height-Angle Charts." Naval
    %       Research Laboratory, NRL Report 6650, Jan. 22, 1968.
    %   [2] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [3] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off., 1959.
    '''
    userSetMaxLims = False
    if (rmax is not None) and (hmax is not None):
        userSetMaxLims = True
        
    rngunit, rngunitDisp, rngfac = unitnamemod(rngunit)
    htunit, htunitDisp, htfac = unitnamemod(htunit)
    addDataTips = False
    
    hVCP = rhaplot(userSetMaxLims, vcp, vcpangles, rngunit, rngunitDisp, rngfac, htunit, htunitDisp, htfac, spow, Ns, expo, anht, fc, ec, hAxes, addDataTips, rmax, hmax)

def unitsratio(to_in, from_in):
    '''
    %UNITSRATIO Unit conversion factors
    %
    %   RATIO = UNITSRATIO(TO, FROM) returns the number of TO units per one
    %   FROM unit.  For example, UNITSRATIO('cm', 'm') returns 100 because
    %   there are 100 centimeters per meter.  UNITSRATIO makes it easy to_in
    %   convert from_in one system of units to_in another.  Specifically, if X is
    %   in units FROM and Y is in units TO, then the conversion equation is
    %
    %                  Y = UNITSRATIO(TO, FROM) * X.
    %
    %   TO and FROM may be any of the length units supported by
    %   validateLengthUnit, or may be one of the following angle units:
    %
    %     Angle Unit            Valid Inputs
    %     ----------            ------------
    %     radian               'rad', 'radian', 'radians'
    %     degree               'deg', 'degree', 'degrees'
    %
    %   Examples
    %   --------
    %   % Approximate mean earth radius in meters
    %   radiusInMeters = earthRadius('meters')
    %   % Conversion factor
    %   feetPerMeter = unitsratio('feet', 'meter')
    %   % Radius in (international) feet:
    %   radiusInFeet = feetPerMeter * radiusInMeters
    %
    %   % The following prints a true statement for any valid TO, FROM pair:
    %   to_in   = 'feet';
    %   from_in = 'mile';
    %   sprintf('There are %g %s per %s.', unitsratio(to_in,from_in), to_in, from_in)
    %
    %   % The following prints a true statement for any valid TO, FROM pair:
    %   to_in   = 'degrees';
    %   from_in = 'radian';
    %   sprintf('One %s is %g %s.', from_in, unitsratio(to_in,from_in), to_in)
    %
    %   See also validateLengthUnit.
    
    % Copyright 2002-2017 The MathWorks, Inc.
    '''
    # Validate units and convert to standard names.
    degreeStrings = ['deg', 'degree', 'degrees']
    radianStrings = ['rad', 'radian', 'radians']
    
    try:
        to_in = validateLengthUnit(to_in, 'UNITSRATIO', 'TO', 1)
    except Exception as e:
        if any(strcmpi(to_in, degreeStrings)):
            to_in = 'degree'
        elif any(strcmpi(to_in, radianStrings)):
            to_in = 'radian'
        else:
            raise e
    
    try:
        from_in = validateLengthUnit(from_in, 'UNITSRATIO', 'FROM', 2)
    except:
        if any(strcmpi(from_in, degreeStrings)):
            from_in = 'degree'
        elif any(strcmpi(from_in, radianStrings)):
            from_in = 'radian'
        else:
            exception.throw()
    
    # Define a relationship graph by specifying a scaling factor for each pair
    # of adjacent units (using the standard names). Each unit on the left is
    # defined in terms of the unit on the right via the factor provided, which
    # is the value that will be returned by unitsratio(RIGHT, LEFT). For
    # example, unitsratio('meter','micron') returns 1e-6.
    graph = [
            ['micron',             'meter',   1e-6],
            ['millimeter',         'meter',   1e-3],
            ['centimeter',         'meter',   1e-2],
            ['kilometer',          'meter',   1e+3],
            ['nautical mile',      'meter',   1852],
            ['foot',               'meter',   0.3048],
            ['mile',               'foot',    5280],
            ['inch',               'foot',    1/12],
            ['yard',               'foot',    3],
            ['U.S. survey foot',   'meter',   1200/3937],
            ['U.S. survey mile',   'U.S. survey foot', 5280],
            ['Clarke''s foot',     'meter',   0.3047972654],
            ['German legal metre', 'meter',   1.0000135965],
            ['Indian foot',        'meter',   12/39.370142],
            ['degree',             'radian',  pi/180],
            ]
    graph = pd.DataFrame(graph)
    
    # Do a depth-first search of the directed graph corresponding
    # to_in the definitions array, recursively searching for a
    # path from_in FROM to_in TO.
    ratio, _ = searchgraph(to_in, from_in, graph, [])

    # A return value of NaN in RATIO indicates that no connection exists.
    if isnan(ratio):
        raise ValueError(f"unitsratio: unable to_in convert units to_in {to_in} from_in {from_in}")

    return unlevel(ratio)
        
def searchgraph(to_in, from_in, graph, history):
    
    # Assume a dead-end unless/until a path is found from_in FROM to_in TO.
    ratio = nan
    
    # Stop here if FROM has already been checked (avoid loops in the graph).
    if np.any(strcmp(from_in, history)):
        return ratio, history
    
    # Append FROM to_in the list of nodes that have been visited.
    if numel(history) > 0:
        history[-1] = from_in
    
    # Find occurrences of FROM and TO in columns 1 and 2 of GRAPH.
    from1 = find(strcmp(from_in, graph.iloc[:, 0].values))
    from2 = find(strcmp(from_in, graph.iloc[:, 1].values))
    to1   = find(strcmp(to_in,   graph.iloc[:, 0].values))
    to2   = find(strcmp(to_in,   graph.iloc[:, 1].values))
    
    # See if there's a direct conversion from_in TO to_in FROM:
    # If there's a row with TO in column 1 and FROM in column 2, then
    # column 3 of that row contains the conversion factor
    # from_in FROM to_in TO.
    i = intersect(to1, from2)
    if numel(i) == 1:
        ratio = 1 / graph.iloc[i, 2].values
        return ratio, history
    
    # See if there's a direct conversion from_in FROM to_in TO:
    # If there's a row with FROM in column 1 and TO in column 2,
    # then column 3 of that row contains the conversion factor.
    i = intersect(to2, from1)
    if numel(i) == 1:
        ratio = graph.iloc[i, 2].values
        return ratio, history
    
    # Recursively search for conversion to_in TO from_in each node adjacent
    # to_in FROM.
    
    # Search from_in the adjacent nodes with a direct conversion _from_
    # FROM.  If a conversion factor (non-NaN) to_in TO is found from_in
    # one of these adjacent nodes, then multiply it by the conversion
    # factor from_in FROM to_in that neighbor (divide by the defining
    # factor in column 3 of GRAPH).
    for i in range(0, numel(from2)):
       n = from2(i)
       ratio, history = searchgraph(to_in, graph.iloc[n, 0].values, graph, history);
       if ~isnan(ratio):
           ratio = ratio / grap.iloc[n, 2].values
           return ratio, history
    
    # Search from_in the adjacent nodes with a direct conversion _to_ FROM.
    # If a conversion factor (non-NaN) to_in TO is found from_in one of these
    # adjacent nodes, then divide it by the conversion factor from_in FROM
    # to_in that neighbor (multiply by the defining factor in column 3).
    for i in range(0, numel(from1)):
       n = from1[i]
       ratio, history = searchgraph(to_in, graph.iloc[n, 1], graph, history)
       if ~isnan(ratio):
           ratio = ratio * graph.iloc[n, 2]
           return ratio, history
    
    return ratio, history

def validateLengthUnit(unit, funcName = None, varName = None, argIndex = None):
    '''
    %validateLengthUnit Validate and standardize length unit name
    %
    %   standardName = validateLengthUnit(UNIT, FUNCNAME, VARNAME, ARGINDEX)
    %   checks that UNIT is a valid length unit and converts it to a standard
    %   unit name.  The following table lists each standard name that is
    %   supported, and the values that correspond to that name. The function is
    %   case-insensitive with respect to its input.  Spaces, periods, and
    %   apostrophes are ignored.  Plural forms are accepted in most cases, but
    %   the result (standardName) is always singular. The optional inputs
    %   FUNCNAME, VARNAME, and ARGINDEX may be included for use in error
    %   message formatting, with behavior identical to that provided by the
    %   VALIDATEATTRIBUTES inputs of the same names.
    %
    %     Standard Name            Supported Names
    %     -------------            ---------------
    %     meter                'm', 'meter(s)', 'metre(s)'
    %
    %     centimeter           'cm', 'centimeter(s)', 'centimetre(s)'
    %
    %     millimeter           'mm', 'millimeter(s)', 'millimetre(s)'
    %
    %     micron               'micron(s)'
    %
    %     kilometer            'km', 'kilometer(s)', 'kilometre(s)'
    %
    %     nautical mile        'nm', 'naut mi', 'nautical mile(s)'
    %
    %     foot                 'ft',   'international ft',
    %                          'foot', 'international foot',
    %                          'feet', 'international feet'
    %
    %     inch                 'in', 'inch', 'inches'
    %
    %     yard                 'yd', 'yds', 'yard(s)'
    %
    %     mile                 'mi', 'mile(s)', 'international mile(s)'
    %
    %     U.S. survey foot     'sf',
    %                          'survey ft',   'US survey ft', 'U.S. survey ft',
    %                          'survey foot', 'US survey foot', 'U.S. survey foot',
    %                          'survey feet', 'US survey feet', 'U.S. survey feet',
    %
    %     U.S. survey mile     'sm', 'survey mile(s)', 'statute mile(s)',
    %     (statute mile)       'US survey mile(s)', 'U.S. survey mile(s)'
    %
    %     Clarke's foot        'Clarke''s foot', 'Clarkes foot'
    %
    %     German legal metre   'German legal metre', 'German legal meter'
    %
    %     Indian foot          'Indian foot'
    %
    %   Examples
    %   --------
    %   % The following return 'foot'
    %   validateLengthUnit('foot')
    %   validateLengthUnit('feet')
    %   validateLengthUnit('international feet')
    %
    %   % The following return 'kilometer'
    %   validateLengthUnit('kilometer')
    %   validateLengthUnit('km')
    %   validateLengthUnit('kilometre')
    %   validateLengthUnit('kilometers')
    %   validateLengthUnit('kilometres')
    %
    %   % The following return 'U.S. survey foot'
    %   validateLengthUnit('U.S. survey foot')
    %   validateLengthUnit('US survey ft')
    %   validateLengthUnit('sf')
    %   validateLengthUnit('U. S. survey feet')
    %   validateLengthUnit('u s survey foot')
    %
    %   % In the following, a non-char input to validateLengthUnit results in
    %   % an error message referencing a function name, 'FOO', a variable name,
    %   % 'UNIT', and an argument number, 5.
    %   validateLengthUnit(17,'FOO','UNIT',5)
    %
    %   See also UNITSRATIO.
    
    % Copyright 2011-2017 The MathWorks, Inc.
    '''
    # Save a copy of the original UNIT, to be used in case of error.
    originalString = unit
    
    # Preprocess the UNIT:
    #   1. Remove periods and apostrophes.
    #   2. Remove spaces.
    #   3. Convert to lower case.
    unit = unit.replace(".", "")
    unit = unit.replace("'", "")
    unit = unit.replace(" ", "")
    unit = unit.lower()
    
    # Synonyms: Values of UNIT found in the first column corresponded to
    # standard names in the second column.  However, some names in the
    # second column (e.g., 'ussurveyfoot') are not fully in standard form,
    # because they have been contracted and converted to lower case (to match
    # the processing applied to UNIT).  Such names are fully standardized in
    # a subsequent filtering step.
    synonyms = np.array([
        [                 'm',  'meter'],
        [            'meters',  'meter'],
        [             'metre',  'meter'],
        [            'metres',  'meter'],
        [           'microns',  'micron'],
        [                'cm',  'centimeter'],
        [       'centimeters',  'centimeter'],
        [        'centimetre',  'centimeter'],
        [       'centimetres',  'centimeter'],
        [                'mm',  'millimeter'],
        [       'millimeters',  'millimeter'],
        [        'millimetre',  'millimeter'],
        [       'millimetres',  'millimeter'],
        [                'km',  'kilometer'],
        [        'kilometers',  'kilometer'],
        [         'kilometre',  'kilometer'],
        [        'kilometres',  'kilometer'],
        [                'nm',  'nauticalmile'],
        [            'nautmi',  'nauticalmile'],
        [     'nauticalmiles',  'nauticalmile'],
        [                'in',  'inch'],
        [            'inches',  'inch'],
        [                'yd',  'yard'],
        [               'yds',  'yard'],
        [             'yards',  'yard'],
        [                'ft',  'foot'],
        [   'internationalft',  'foot'],
        [ 'internationalfoot',  'foot'],
        [              'feet',  'foot'],
        [ 'internationalfeet',  'foot'],
        [                'mi',  'mile'],
        [             'miles',  'mile'],
        [ 'internationalmile',  'mile'],
        ['internationalmiles',  'mile'],
        [                'sf',  'ussurveyfoot'],
        [        'surveyfoot',  'ussurveyfoot'],
        [          'surveyft',  'ussurveyfoot'],
        [        'ussurveyft',  'ussurveyfoot'],
        [      'ussurveyfoot',  'ussurveyfoot'],
        [        'surveyfeet',  'ussurveyfoot'],
        [      'ussurveyfeet',  'ussurveyfoot'],
        [                'sm',  'ussurveymile'],
        [        'surveymile',  'ussurveymile'],
        [       'surveymiles',  'ussurveymile'],
        [       'statutemile',  'ussurveymile'],
        [      'statutemiles',  'ussurveymile'],
        [     'ussurveymiles',  'ussurveymile'],
        [       'clarkesfoot',  'clarkesfoot'],
        [  'germanlegalmeter',  'germanlegalmetre'],
        [        'indianfoot',  'indianfoot'],
    ])

    if ~np.any(strcmp(unit, synonyms[:, 1])):
        # UNIT doesn't match any of the standard length unit names in the
        # second column of the synonyms array; see if it matches any of the
        # synonyms in the first column. If a match is found, replace the
        # synonym with the standard name.
        index = strcmp(unit, synonyms[:, 0]);
        if any(index):
            # The first column is known to be free of duplicates, so we
            # can assume that exactly one element in index is true.
            # Select the corresponding unit from the second column.
            unit = synonyms[index, 1]
        else:
            raise ValueError(f"validateLengthUnit: unsupported unit: {originalString}")
    
    # Filter contracted names: UNIT now corresponds to a standard name, except
    # that missing spaces, missing periods, and capitalization may need to be
    # restored.
    contractedNames = np.array([
        ['nauticalmile',     'nautical mile'],
        ['ussurveyfoot',     'U.S. survey foot'],
        ['ussurveymile',     'U.S. survey mile'],
        ['clarkesfoot',      'Clarke''s foot'],
        ['germanlegalmetre', 'German legal metre'],
        ['indianfoot',       'Indian foot'],
    ])
    index = strcmp(unit, contractedNames[:, 0])
    if np.any(index):
        # There's a (unique) match with an element of the first column of
        # contractedNames; replace unit with the corresponding element from the
        # second column (the uncontracted form).
        unit = contractedNames[index, 1]
    
    return unit

def unitnamemod(unitNameDisp):
    '''
    %unitnamemod Updates the unit name for the units ratio function 
    %   [UNITNAME,FAC] = unitnamemod(UNITNAME) converts UNITNAME to the
    %   appropriate unit name for use by unitratio, as well as an additional
    %   conversion factor FAC in the case of kft. 
    %   
    %   Unit name changes:
    %      Input Name     |     Output Name     |     Conversion Factor F
    %      kft                  ft                    1e3
    %      nmi                  nm                    1
    %
    %   Allowed units are the same as unitsratio with the added kft option.
    '''
    unitNameToConvert = ['nmi', 'kft']
    unitNameConversion = ['nm', 'ft']
    facConversion = [1, 1e3]
    
    tf = [x.lower() == unitNameDisp.lower() for x in unitNameToConvert]
    
    if any(tf):
        unitNameOut = unitNameConversion[tf]
        fac = facConversion[tf]
    else:
        unitNameOut = unitNameDisp
        fac = 1

    return unitNameOut, unitNameDisp, fac
    
def rhaplot(
    userSetMaxLims,
    vcp,
    vcpangles,
    rngunit,
    rngunitDisp,
    rngfac,
    htunit,
    htunitDisp,
    htfac,
    spow,
    Ns,
    expo,
    anht,
    fc,
    ec,
    hAxes,
    addDataTips,
    rmax,
    hmax,
):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rhaplot Range-height-angle plot (Blake chart)
    %   rhaplot(userSetMaxLims,vcp,rngunit,rngunitDisp,rngfac,htunit, ...
    %   htunitDisp,htfac,spow,Ns,expo,anht,fc,ec,hAxes,addDataTips)
    %
    %   userSetMaxLims = Logical indicating whether user set rmax and hmax
    %   vcp            = Vertical coverage contour (in rngunit)
    %   vcpangles      = Vertical coverage angles (in degrees)
    %   rngunit        = Range unit of x-axis and vcp used for conversions
    %   rngunitDisp    = Range unit to display in plot
    %   rngfac         = Additional conversion factor to account for kft
    %   htunit         = Height unit of y-axis and anht used for conversions
    %   htunitDisp     = Height unit to display in plot
    %   htfac          = Additional conversion factor to account for kft
    %   spow           = Scale power between 0 and 1
    %   Ns             = Surface refractivity in N-units
    %   expo           = Refraction exponent in 1/km
    %   anht           = Antenna height (in htunit)
    %   fc             = Face (patch) color, N x 3, where N can be 1 or the
    %                    number of columns in vcp
    %   ec             = Edge (contour) color, N x 3, where N can be 1 or the
    %                    number of columns in vcp
    %   hAxes          = Plot axes
    %   addDataTips    = Logical indicating whether data tips should be added.
    %                    Turn off to increase speed.
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also blakechart, radarDesigner.
    
    %   Copyright 2021-2024 The MathWorks, Inc.
    '''

    # Nested helper functions
    # --------------------------
    def getheightXY():
        hts = get_heightsamples(hmax, E, anht) * unitsratio('km', htunit)
        srange = slantrange(hts, theta_0, Ns, expo, anhtkm)
        htX, htY = rngth2xy(srange, theta_1, rmaxkm, spow, E)
        
        return htX, htY, hts
    
    def getboundingxy(bdhtx, bdhty, srGeRmax):
        # Solves for the xy coordinate of the boundary of the chart
        if srGeRmax:
            # When range corresponding to hmax exceeds rmax
            # Solve for point where loci for const hmax intersects with const rmax
            if (hmaxkm == rmaxkm): # fzero fails to solve in this case
                intrth0 = pi / 2
            else: # Use fzero to solve
                # make sure answer is between 0 and pi/2
                thsol = brentq(lambda th: slantrange(hmaxkm, th, Ns, expo, anhtkm) - rmaxkm, 0, pi / 2)
                intrth0 = rem(abs(thsol), pi/2)

            intrth1 = atan(E * tan(intrth0)) # Angle on chart
            _, yintr = rngth2xy(rmaxkm, intrth1, rmaxkm, spow, E)
            bdx = [cos(linspace(0, intrth0)), bdhtx[bdhty > yintr],  0]
            bdy = [E * sin(linspace(0, intrth0)), bdhty[bdhty > yintr], 0]
        else:
            # Solve for angle where the const range curve meets y = 1 line
            intrth0 = brentq(lambda th: (E * sin(th) - 1), 0, pi / 2)
            crngth = [theta_0[theta_0 < intrth0], intrth0]
            bdx = [cos(crngth), 0, 0]
            bdy = [E * sin(crngth), 1, 0]

        return bdx, bdy, intrth0
    
    def getangleXY(srGeRmax):
        # Returns X,Y coordinates for angle spokes
        # Set index for angle to be marked
        # theta_0 is 181 samples
        if E <= 5:
            idxa = np.hstack([1, list(range(101, 171, 10)), 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 5 and E <= 25):
            idxa = np.hstack([list(range(1, 101, 10)), list(range(111, 151, 10)), 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 25 and E <= 50):
            idxa = np.hstack([list(range(1, 101, 10)), 111, 121, 136, 151, 171, 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 50 and E<= 75):
            idxa = np.hstack([list(range(1, 101, 10)), 111, 121, 136, 151, 181]) - 1 # may need to adjust these for 0 indexing
        else:
            idxa = np.hstack([list(range(1, 11)), list(range(21, 101, 10)), 111, 181]) - 1 # may need to adjust these for 0 indexing
        
        if srGeRmax:
            theta_0s = theta_0[idxa]
            tha0 = theta_0s[theta_0s <= intrth0] # part with const rng boundary
            st = find(theta_0s > intrth0, 1, 'first')
            thb0 = theta_0[idxa[unlevel(st):]]
            thsz = length(tha0) + length(thb0)
            
            angX = np.vstack([[zeros(thsz)], np.hstack([cos(tha0), htX[-1, idxa[unlevel(st):]]])])
            angY = np.vstack([[zeros(thsz)], np.hstack([E * sin(tha0), htY[-1, idxa[unlevel(st):]]])])
        else:
            theta_0s = theta_0[idxa] # Sampled theta_0
            tha0 = theta_0s[theta_0s <= intrth0] # part with const rng boundary
            thb0 = theta_0s[theta_0s > intrth0] # part with y = 1 boundary
            thb1 = atan(E * tan(thb0)) # Convert to angles on chart for this portion
            thsz = length(tha0) + length(thb1)
            
            angX = np.vstack([[zeros(thsz)], np.hstack([cos(tha0), 1 / tan(thb1)])])
            angY = np.vstack([[zeros(thsz)], np.hstack([E * sin(tha0), ones(1, length(thb1))])])

        return angX, angY
    
    def getrangeXY(srGeRmax):
        # Range related calculations
        rngs = get_rangesamples(rmax, E) * unitsratio('km', rngunit)
        rngX = cos(singcolvec(theta_0)) * ((rngs/rmaxkm) ** spow)
        rngY =  E * sin(singcolvec(theta_0)) * ((rngs / rmaxkm) ** spow)
        
        if srGeRmax:
            # Const range
            for col in range(0, len(rngs)):
                chtbdY = interp1(htX[-1, :], htY[-1, :], rngX[:, col],'pchip') # Boundary # -1's were end
                temp = rngY[:, col]
                temp[temp > chtbdY] = chtbdY[temp > chtbdY]
                rngY[:, col] = temp
        else:
            # Const Range
            rngY[rngY > 1] = 1

        return rngX, rngY, rngs
    
    curveColor = '--mw-graphics-borderColor-axes-tertiary'
    
    # Determine max range and height
    rmax, rmaxkm, rmaxPlotkm, hmax, hmaxkm, hmaxPlotkm = rmaxhmax(userSetMaxLims, vcp, vcpangles, rngunit, rngfac, htunit, htfac, Ns, expo, anht, rmax, hmax)

    # Convert to SI units. All computation done in SI units
    vcpkm = singcolvec(vcp * unitsratio('km', rngunit))
    vcpangles = deg2rad(vcpangles)
    anhtkm = anht * unitsratio('km', htunit)
    
    # Compute ellipticity factor
    E = (rmaxkm / hmaxkm) ** spow # Ref 1, eqn 8
    
    # Define el-angle range (true angle)
    theta_0 = get_elangle() # True angle
    theta_1 = true2chart(theta_0, E) # Angle on chart
    
    # Computation begins
    srGeRmax = slantrange(hmaxkm, 0, Ns, expo, anhtkm) >= rmaxkm
    htX, htY, hts = getheightXY()
    bdx, bdy, intrth0 = getboundingxy(htX[-1, :], htY[-1, :], srGeRmax) # -1's were end
    rngX, rngY, rngs = getrangeXY(srGeRmax)
    angX, angY = getangleXY(srGeRmax)
    
    # Set axes
    xlimMax = (rmaxPlotkm / rmaxkm) ** spow
    ylimMax = (hmaxPlotkm / hmaxkm) ** spow

    fig, hAxes = plt.subplots(figsize=(8, 5))
    # plt.figure(figsize=(8, 5))
    plt.axis([0, xlimMax, 0, ylimMax])
    # set(hAxes, 'Color', 'none')
    # hold(hAxes, 'on')
    
    # # Plot background patch
    # patchBackgroundColor = '--mw-graphics-backgroundColor-axes-primary'
    # patchBorderColor = '--mw-graphics-borderColor-axes-primary'
    # hBackground = patch(hAxes, bdx, bdy, 'w', 'HitTest', 'off', 'Tag', 'background')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hBackground, 'FaceColor', patchBackgroundColor)
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hBackground, 'EdgeColor', patchBorderColor)
    # hBackground.Annotation.LegendInformation.IconDisplayStyle = 'off'
    
    # Find constant height x, y points
    if srGeRmax:
        htXGtCosTheta0 = np.array(bsxfun("gt", htX, cos(theta_0))).astype(int)
        idxFind = np.pad(singrowvec(matlabsum(htXGtCosTheta0, 2) >= 1), (0, numel(htY) - numel(singrowvec(matlabsum(htXGtCosTheta0, 2) >= 1))))
        rowVec = np.array(list(range(0, numel(htY))))
        
        for row in rowVec[idxFind]: # Each row -> each const ht curve
            # For each height curve
            idxc = find(htXGtCosTheta0[row, :], 1, 'last')
            
            # Interpolation to make the plot smooth at edge
            htX[row, idxc] = cos(theta_0[idxc])
            htY[row, idxc] =  interp1(htX[row, :], htY[row, :], cos(theta_0[idxc]), 'pchip')
            htY[row, :int(idxc - 1)] = nan
    
    # # Plot custom range, height, angle axes
    plotCustomAxes(hAxes, rngX, rngY, htX, htY, angX, angY, curveColor)
    
    # Create angle labels
    tx = angX[1, :]
    ty = angY[1, :]
    angth, angr = cart2pol(tx, ty)
    angr = angr + 0.01 # Position the label slightly outside the boundary
    angx, angy = pol2cart(angth, angr)
    
    # Add angle labels at edge of plot
    anglbl = rad2deg(chart2true(angth, E))
    formats = ['{:.0f}', '{:.1f}']
    numAng = len(angth)
    numFormat = [1] * numAng
    numFormat = [2 if (anglbl[i] < 1 and anglbl[i] != 0) else 1 for i in range(numAng)]
    angStr = [formats[numFormat[i]-1].format(anglbl[i]) + chr(176) for i in range(numAng)]
    hAngText = text(hAxes, angx, angy, angStr, VerticalAlignment = 'baseline', Clipping = 'on')
    # fixOverlappingText(hAngText)
    
    # Add angle labels within the boundary
    if (xlimMax < 1) or (ylimMax < 1): # #ok<BDSCI>
        maxR = sqrt(xlimMax ** 2 + ylimMax ** 2)
        angr = maxR - 0.3 # Position label slightly inside the boundary
        angx, angy = pol2cart(angth, angr)
        hAngInnerText = text(hAxes, angx, angy, angStr, VerticalAlignment = 'baseline', Clipping = 'on', BackgroundColor = 'none', Tag = 'AngleLabels')
        # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngInnerText,'Color',curveColor)
        # fixOverlappingText(hAngInnerText)
    
    # Translate range/angle to Cartesian on chart
    _, c = size(vcpkm)
    vcpX = Nan(numel(vcpangles), c)
    vcpY = Nan(numel(vcpangles), c)
    vcpangc = true2chart(vcpangles, E)
    idxPlot = vcpangc >= 0
    for idx in range(c - 1, -1, -1):
        vcpX[idxPlot, idx], vcpY[idxPlot, idx] = rngth2xy(vcpkm[idxPlot, idx], vcpangc[idxPlot], rmaxkm, spow, E)
    
    # Plot pattern patch
    for idx in range(c - 1, -1, -1):
        idxPlot = ~isnan(vcpX[:, idx]) & ~isnan(vcpY[:, idx])
        # hPatch = patch(hAxes, [[0], [vcpX[idxPlot, idx]], [0]], [[0], [vcpY[idxPlot, idx]], [0]], [1, 1, 1], 'EdgeColor', 'none', 'LineStyle', 'none', 'HitTest', 'off', 'FaceAlpha', 0.6, 'Tag', sprintf('vcpPatch-%d', idx))
        # applyPatchColor(hPatch, fc, idx)
    
    # Plot contour lines
    h_vcp = np.empty(c, dtype=object) 
    for idx in range(c - 1, -1, -1):
        # set endpoints to zero to close fill
        vcpY[:, idx][0] = 0
        vcpY[:, idx][-1] = 0
        vcpX[:, idx][0] = 0
        vcpX[:, idx][-1] = 0
        
        h_vcp[idx] = plt.plot(vcpX[:, idx], vcpY[:, idx])#, 'Tag', sprintf('vcpLine-%d', idx))
        hAxes.fill(vcpX[:, idx], vcpY[:, idx], color='C0', alpha=0.5)

        # applyPlotColor(h_vcp[idx], ec, idx)
    
        # # Plot pattern with no interactions. Interactions will be on custom
        # # data tips.
        # h_vcp(idx).PickableParts = 'none'
    
    #----------- X and Y-Axis Ticks and Labels -----------------------------------
    # Range and height axis labels
    xnorm = (rngs / rmaxkm) ** spow
    ynorm = (hts / hmaxkm) ** spow
    if E <= 5:
        tidx = np.hstack([1, 15, list(range(19, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([1, 20, 25, list(range(29, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing
    elif (E > 5 and E < 50):
        tidx = np.hstack([1, list(range(10, 18, 2)), list(range(19, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([1, list(range(10, 18, 2)), list(range(19, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing
    else:
        tidx = np.hstack([list(range(1, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([list(range(1, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing

    tidx = tidx - 1
    tidy = tidy - 1

    tidx, _ = unique(tidx)
    tidy, _ = unique(tidy)
    
    XTicks = xnorm[tidx]
    XTickLabel = rngs[tidx] * (1 / rngfac) * unitsratio(rngunit, 'km')
    # XTicks, XTickLabel = fixOverlappingTicks(hAxes, hAngText, XTicks, XTickLabel, 'x')
    
    YTicks = ynorm[tidy]
    YTickLabel = hts[tidy] * (1 / htfac) * unitsratio(htunit, 'km')
    # YTicks, YTickLabel = fixOverlappingTicks(hAxes, hAngText, YTicks, YTickLabel, 'y')
    
    # # Set X tick marks
    # set(hAxes, 'XTick', XTicks)
    # set(hAxes, 'XTickLabel', XTickLabel)
    hAxes.set_xlabel(f"Range ({rngunitDisp})")
    
    # # Set Y tick marks
    # set(hAxes, 'YTick', YTicks)
    # set(hAxes, 'YTickLabel', YTickLabel)
    hAxes.set_ylabel(f"Height ({htunitDisp})")
    
    # # Add title
    title_text_object = hAxes.set_title("Blakechart")
    pos = title_text_object.get_position()
    hAxes.set_title("Blakechart", x = pos[0], y = pos[1] + 0.04)
    
    return h_vcp

# Validation/parsing helper functions
# -----------------------------------------+
def rmaxhmax(
    userSetMaxLims,
    vcp,
    vcpangles,
    rngunit,
    rngfac,
    htunit,
    htfac,
    Ns,
    expo,
    anht,
    rmax,
    hmax,
):
    # Set rmax/hmax
    if userSetMaxLims:
        # Use set max lims
        validateboundarymax(rmax, hmax)
        
        # Adjust for potential kft input
        rmaxPlot = rmax * rngfac
        hmaxPlot = hmax * htfac
        
        # Estimate actual maximum range and height
        rmax, hmax = estrmaxhmax(vcp, vcpangles, anht, rngunit, htunit, Ns, expo, rmax, hmax)
    else:
        # Estimate maximum range and height
        rmax, hmax = estrmaxhmax(vcp, vcpangles, anht, rngunit, htunit, Ns, expo)
        rmaxPlot = rmax
        hmaxPlot = hmax
    
    # Convert to km
    rmaxkm = rmax * unitsratio('km', rngunit)
    hmaxkm = hmax * unitsratio('km', htunit)
    rmaxPlotkm = rmaxPlot * unitsratio('km', rngunit)
    hmaxPlotkm = hmaxPlot * unitsratio('km', htunit)
    anhtkm = anht * unitsratio('km',htunit)
    hmaxkm = hmaxkm + anhtkm
    hmaxPlotkm = hmaxPlotkm + anhtkm
    
    # RMAX >= HMAX
    if rmaxkm < hmaxkm:
        rmaxkm = hmaxkm
        
    return rmax, rmaxkm, rmaxPlotkm, hmax, hmaxkm, hmaxPlotkm

def validateboundarymax(rmax, hmax):
    validateattributes(rmax,{'double'},{'positive','nonzero','scalar','finite', 'real'}, 'blakechart', 'RMAX')
    validateattributes(hmax,{'double'},{'positive','nonzero','scalar','finite', 'real'}, 'blakechart', 'HMAX')
    
# Calculation helper functions
# -----------------------------------------+
def slantrange(height, theta, Ns, d, anht):
    '''
    %SLANTRANGE Coordinate on range-height-angle chart
    %   SRANGE = SLANTRANGE(HEIGHT, THETA) Returns propagated range (SRANGE)
    %   for given HEIGHT and ANGLE on a range-height-angle chart
    %   HEIGHT Height from earth surface in km
    %
    %   Ref:
    %   [1] Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
    %       1968 Eq 4, p 11.
    %   [2] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968.
    %   [3] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''
    
    # All calculation are done in SI units
    if (numel(height) == 1) and (numel(theta) == 1):
        srange = rangeIntegralCRPL(height, anht, theta, Ns, d)
    else:
        hlog = logspace(-12, 2, int(1e4))
        hint = matlabunion(height, hlog)
        y = rangeIntegrandCRPL(hint, anht, theta, Ns, d)
        srangeCumtrapz = cumtrapz(hint, y)
        idxHeight = [int(i) for i in interp1(hint, list(range(0, numel(hint))), height, 'nearest')]
        srange = srangeCumtrapz[idxHeight, :]
        srange[:, 0] = rangeIntegralCRPL(height, anht, theta[0], Ns, d)

    return srange

def rangeIntegrandCRPL(h, ha, theta0, Ns, k):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rangeIntegrandCRPL Integrand for propagated range using CRPL 
    %   R = rangeIntegrandCRPL(H,HA,THETA0,NS,K) calculates the propagated
    %   range R in km using the CRPL atmospheric model. 
    %
    %   The inputs are as follows: 
    %      H      = Target height (km)
    %      HA     = Antenna height (km)
    %      THETA0 = Elevation angle at radar (rad)
    %      NS     = Surface refractivity (N-units, 0 height)
    %      K      = CRPL decay constant (1/km)
    %      
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also height2range, height2grndrange, el2height.
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''

    h = np.array(singcolvec(h))
    theta0 = np.array(singrowvec(theta0))
    
    r0 = Rearth * 1e-3 + ha # Eqn 8 (km)
    Ns = Ns * 1e-6
    rho0 = Ns * exp(-k * ha) # Eqn 7
    n = 1 + rho0 * exp(-k * h)
    u = (1 + rho0) ** 2 * sin(theta0) ** 2 - 2 * rho0 - rho0 ** 2
    v = 2 * rho0 * exp(-k * h) + rho0 ** 2 * exp(-2 * k * h)
    w = 2 * h / r0 + h ** 2 / r0 ** 2

    R = np.absolute(n ** 2 * (1 + h / r0) / sqrt(u + v + w + v * w)) # Eqn 16
    
    return R

def rangeIntegralCRPL(heightsIn, haIn, theta0In, Ns, k):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rangeIntegrandCRPL Integrand for propagated range using CRPL 
    %   R = rangeIntegrandCRPL(H,HA,THETA0,NS,K) calculates the propagated
    %   range R in km using the CRPL atmospheric model. 
    %
    %   The inputs are as follows: 
    %      H      = Target height (km)
    %      HA     = Antenna height (km)
    %      THETA0 = Elevation angle at radar (rad)
    %      NS     = Surface refractivity (N-units, 0 height)
    %      K      = CRPL decay constant (1/km)
    %      
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also height2range, height2grndrange, el2height.
    
    %   Copyright 2021-2022 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''
    
    # Expand vectors as necessary
    maxNum, _ = matlabmax([numel(heightsIn), numel(haIn), numel(theta0In)])
    heights = expandVector(heightsIn, maxNum)
    ha = expandVector(haIn, maxNum)
    theta0 = expandVector(theta0In, maxNum)

    if maxNum == 1:
        heights = np.array([heights])
        ha = np.array([ha])
        theta0 = np.array([theta0])
    
    # Only perform integral for values of theta > 0.1 degrees and heights > 1
    # ft (3.048e-04 km). See pg 89 in reference 1.
    idx = np.array(list(range(0, maxNum)))
    idxNear0 = (abs(theta0) < deg2rad(0.1))
    idxIntegrate = ~(idxNear0 & (abs(heights) < 3.048e-04))
    R = zeros(maxNum)
    if np.any(idxIntegrate):
        # Loosen integration error bounds in the case when angle is near 0
        if np.any(idxNear0):
            thisIdx = idxIntegrate & idxNear0
            R[thisIdx] = np.absolute(arrayfun(lambda x: quadgk(lambda h: rangeIntegrandCRPL(h, ha[x], theta0[x], Ns, k), 0, heights[x], AbsTol = 1e-2, RelTol = 0)[0], idx[thisIdx])) # Eqn 2
        
        # Perform normal calculations for larger angles
        if np.any(~idxNear0):
            thisIdx = idxIntegrate & ~idxNear0
            R[thisIdx] = np.absolute(arrayfun(lambda x: quadgk(lambda h: rangeIntegrandCRPL(h, ha[x], theta0[x], Ns, k), 0, heights[x])[0], idx[thisIdx])) # Eqn 2
    
    # Perform approximation outlined in section 5 of reference 1 for values
    # near theta0 = 0 and h = 0
    if np.any(~idxIntegrate):
        R[idx[~idxIntegrate]] = rangeIntegralApprox(heights[~idxIntegrate], ha[~idxIntegrate], theta0[~idxIntegrate], Ns, k) # Eqn 22

    return R
    
def expandVector(xin, maxNum):
    if isscalar(xin):
        xout = xin * ones((1, maxNum))
        # xout = xin(1).*ones(1,maxNum);
    else:
        xout = xin

    return unlevel(xout)
    
def rangeIntegralApprox(h, ha, theta0, Ns, k):
    r0 = Rearth * 1e-3 + ha
    rho0 = Ns * 1e-6 * exp(-k * ha)
    gammaConst = rho0 * k / (1 + rho0) # Eqn 23
    A_R = 1 + rho0 # Eqn 24
    a_R, _ = matlabmax(k * gammaConst - 7 * gammaConst ** 2 + 8 * gammaConst / r0 - 3 / r0 ** 2 + sin(theta0) ** 2 * (10 * gammaConst ** 2 - 8 * gammaConst / r0 - 2 * k * gammaConst + 3 / r0 ** 2), eps()) 
    b_R = 2 * ((2 * gammaConst - 1 / r0) * sin(theta0) ** 2 - gammaConst + 1 / r0)
    c_R = sin(theta0) ** 2
    R = taylorExpansionApprox(h, A_R, a_R, b_R, c_R)
    
    return R
    
def taylorExpansionApprox(h1, A, a, b, c):
    h1 = abs(h1)
    x = (2 * a * h1 - 2 * sqrt(a * c) + 2 * sqrt(a * (a * h1 ** 2 + b * h1 + c))) / (b + 2 * sqrt(a * c))
    F = abs(A / sqrt(a) * log(1 + x))

    return F

    
def rngth2xy(rng, th1, rmaxkm, pow_in, E):
    # RNGTH2XY Maps range and angle th1 (angle on chart) to x, y coordinates
    #   Ref: Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
    #   1968 Eq 4, p 11
    parma = (rng / rmaxkm) ** pow_in # Eqn 11, p 13
    parmb = E * parma # Eqn 12, p 13
    
    # Compute x, y coordinates
    xcor = bsxfun("rdivide", parma, sqrt(1 + (tan(th1) / E) ** 2)) # Eqn 9, p 13
    ycor = bsxfun("times", parmb, sqrt(1 - (xcor / parma) ** 2)) # Eqn 10, p 13
    
    return xcor, ycor
    
def true2chart(theta_tr, E):
    # Converts true angles to angles on chart
    # theta_ch = theta_tr if E = 1
    
    # Eq 13, p 13, NRL Report 6650
    
    theta_ch = atan(E * tan(theta_tr))
    return theta_ch
    
def chart2true(theta_ch, E):
    # Converts true angles to angles on chart
    # theta_ch = theta_tr if E = 1

    # Eq 13, p 13, NRL Report 6650
    
    theta_tr = atan(tan(theta_ch) / E)
    return theta_tr
    
# Sample generation helper functions
# -----------------------------------------+
def estrmaxhmax(vcp, vcpang, anht, rngunit, htunit, Ns, expo, rinput = [], hinput = []):
    # Estimate rmax and hmax limits for blakechart based on values of vcp
    # rmax is estimated using the maximum range covered of vcp
    # hmax is estimated using the mean range coverage of vcp
    # This function needs to be optimized to better estimate hmax.

    if np.ndim(vcp) <= 1:
        vcp = np.array([vcp])

    # Set rmax
    vcpm = vcp * unitsratio('m', rngunit)
    rlim1, _ = matlabmax(vcp[-1, :] * cosd(vcpang)) # Estimate 1 of maximum range 
    rlim1 = nextfactor(rlim1)
    _, idxMaxEst = matlabmax(vcp[-1, :] * sind(vcpang)) # Determine the index of the potential maximum height
    _, idxMinEst = matlabmin(vcp[-1, :] * sind(vcpang)) # Determine the index of the potential minimum height

    rEst = [vcpm[-1, i] for i in [idxMaxEst, idxMinEst]]
    angEst = [vcpang[i] for i in [idxMaxEst, idxMinEst]]
    anhts = repmat(anht * unitsratio('m', htunit) * 1e-3, 1, 2) 
    tgthts = (rEst[:] * sind(angEst[:])).T * 1e-3
    tgthts = tgthts - anhts
    rlim2 = abs(rangeIntegralCRPL(tgthts, anhts, deg2rad(angEst).T, Ns, expo)) * 1e3 * unitsratio(rngunit, 'm')
    rlim2 = nextfactor(matlabmax(rlim2[:])[0])

    if rinput != []:
        rmax, _ = matlabmax(np.vstack([[rlim1], [rlim2], [rinput]]))
    else:
        rmax, _ = matlabmax(np.vstack([[rlim1], [rlim2]]))

    if numel(rmax) > 1:
        rmax[isnan(rmax)] = 100e3 * unitsratio(rngunit, 'm')
    elif isnan(rmax):
        rmax[isnan(rmax)] = 100e3 * unitsratio(rngunit, 'm')
    
    # Set hmax
    hEst = range2height(vcpm[-1, idxMaxEst], anht * unitsratio('m', htunit), vcpang[idxMaxEst], method = 'CRPL', Ns = Ns, rexp = expo, maxNumIter = 1, tol = 1e-4) * unitsratio(htunit, 'm')
    hEst = nextfactor(hEst)
    hmax, _ = matlabmax(np.vstack([singcolvec(hEst), singcolvec(hinput)]))
    hmax = np.array(hmax)
    hmax[isnan(hmax)] = 100e3 * unitsratio(htunit, 'm')
    
    return rmax, hmax
    
def nextfactor(x):
    # Returns the next 10 multiple after x
    # Eg nextfactor(13) returns 20, nextfactor(123) returns 200
    ex = fix(log10(x))
    fac = 10 ** ex
    re = rem(x, fac)
    y = (x - re) / fac
    nf = (y + 1) * fac
    return nf
    
def get_elangle():
    # Returns properly sample elevation angles in deg for RHA chart generation
    # Based on code snippet from NRL report 7098
    # Needs optimization
    idx = 0
    elangle = zeros(181)
    stedIdx = np.array([[0, 100], [10, 91]])
    for p in range(0, 2):
        for k in range(stedIdx[p, 0], stedIdx[p, 1]):
            n = k * 10 ** p
            elev = n / 10
            elangle[idx] = deg2rad(elev)
            idx = idx + 1

    return elangle
    
def get_heightsamples(hmax, E, anht):
    # Returns height sample points in range [0 hmax] for constant height curves
    # Sample points for height scale
    expe = fix(log10(hmax))
    
    if E <= 5:
        exps = expe - 3 # 3 decade scale
    elif (E > 5 and E < 50):
        exps = expe - 2
    else:
        exps = expe - 1
    
    ht0 = float(10) ** exps # Starting point
    ht = zeros(9 * (expe-exps) + 1)
    ht[0] = ht0
    for k in range(exps, expe):
        ht[(k - exps) * 9 + np.arange(1, 10)] = np.linspace(2 * 10 ** k, 10 ** (k + 1), 9)

    htx = ht[-1] + 10 ** (k + 1)
    htsamp = unique(np.hstack([ht, np.arange(htx, hmax + 1, 0.5 * 10 ** (k + 1)), hmax]))[0]
    htsamp = htsamp + anht

    return htsamp
    
def get_rangesamples(rmax, E):
    # Returns range sample points in range [0 rmax] for constant height curves
    # Sample points for height scale
    expe = fix(log10(rmax))
    
    if E <= 5:
        exps = expe - 3 # 3 decade scale
    elif (E > 5 and E < 50):
        exps = expe - 2 # 2 decade
    else:
        exps = expe - 1 # 1 decade
    
    rng0 = float(10) ** exps # Starting point
    rng = zeros(9 * (expe - exps) + 1)
    rng[0] = rng0
    for k in range(exps, expe):
        rng[(k - exps) * 9 + np.arange(1, 10)] = np.linspace(2 * 10 ** k, 10 ** (k + 1), 9)

    rngt = rng[-1]  + 0.5 * 10 ** (k + 1)
    rngsamp = unique(np.hstack([rng, np.arange(rngt, rmax + 1, 10 ** (k + 1)), rmax]))[0]

    return rngsamp

def range2height(
    R,
    anht,
    el,
    method = "Curved", # Method
    Re = (4 / 3) * Rearth, # EffectiveEarthRadius
    Ns = 313, # SurfaceRefractivity
    rexp = 0.143859, # RefractionExponent
    maxNumIter = None, # MaxNumIterations
    tol = 1e-6, # Tolerance
):
    '''
    %range2height Convert propagated range to target height
    %   TGTHT = range2height(R,ANHT,EL) returns the target height in meters
    %   assuming a curved Earth model with a 4/3 effective Earth radius, which
    %   is an approximation used for modeling refraction effects in the
    %   troposphere. R is the propagated range between the target and the
    %   sensor in meters. ANHT is the sensor height in meters. EL is the
    %   initial elevation angle of the ray leaving the sensor, also known as
    %   the local elevation angle, in degrees. Heights are referenced to the
    %   ground.
    %
    %   R, ANHT, and EL can either be scalars or matching M-length vectors.
    %
    %   TGTHT = range2height(...,'Method',METHOD) specifies the model for the
    %   calculation as 'Curved' | 'Flat' | 'CRPL'. The default method is
    %   'Curved'.
    %      - 'Curved'
    %         Assumes a curved Earth model with a 4/3 effective Earth radius,
    %         which is an approximation used for modeling refraction effects in
    %         the troposphere. Alternative values for the effective Earth
    %         radius can be set using the 'EffectiveEarthRadius' name-value
    %         argument. This is the default method.
    %      - 'Flat'
    %         Assumes a flat Earth model. In the case of a flat Earth model,
    %         the effective Earth radius is infinite. 
    %      - 'CRPL'
    %         Assumes a curved Earth model where the atmosphere is defined by
    %         the CRPL exponential reference atmosphere with a refractivity of
    %         313 N-units and a refraction exponent (decay constant) of
    %         0.143859 1/km.
    %
    %   TGTHT = range2height(...,'Method','Curved','EffectiveEarthRadius',RE)
    %   specifies the effective Earth radius as a positive scalar in meters.
    %   This input applies only when 'Method' is set to 'Curved' and is
    %   otherwise ignored. The default calculates the effective Earth radius
    %   using a refractivity gradient of -39e-9, which results in approximately
    %   4/3 of the real Earth radius.
    %
    %   TGTHT = range2height(...,'Method','CRPL','SurfaceRefractivity',NS,
    %   'RefractionExponent',REXP) specifies the surface refractivity NS (in
    %   N-Units) as a non-negative scalar and the refraction exponent factor
    %   (decay constant) REXP (in 1/km) as a non-negative scalar for the CRPL
    %   exponential reference atmosphere model. This input is only applicable
    %   when 'Method' is set to 'CRPL' and is otherwise ignored. The default
    %   value of NS is 313 N-Units, and the default value of REXP is 0.143859
    %   1/km.
    %
    %   The atmospheric refraction model is given by
    %        n(h) = 1 + NS*1e-6*exp(-REXP*h),
    %   where h is the height in kilometers and n is the index of refraction.
    %
    %   TGTHT = range2height(...,'MaxNumIterations',MAXITER) specifies the
    %   maximum number of iterations for the CRPL method as a non-negative
    %   scalar integer. This input acts as a safe guard to prevent endless
    %   iterative calculations. This input applies only when 'Method' is set
    %   to 'CRPL' and is otherwise ignored. Defaults to 10.
    %
    %   MAXITER equal to 0 indicates that a quicker but less accurate
    %   non-iterative CRPL calculation should be performed. The non-iterative
    %   calculation was found to have a maximum height error of 0.056388 m
    %   (0.185 ft) at a target height of 30480 m (100,000 ft) and an elevation
    %   angle equal to 0. The height error for the non-iterative method
    %   decreases with decreasing target height and increasing elevation angle.
    %
    %   TGTHT = range2height(...,'Tolerance',TOL) specifies the numerical
    %   tolerance for the CRPL method at which the iterative process is
    %   terminated. This input applies only when 'Method' is set to 'CRPL'
    %   and 'MaxNumIterations' is greater than 0. Otherwise, this input is
    %   ignored. Defaults to 1e-6.
    %
    %   % Examples:
    %
    %   % Example 1: 
    %   %   Determine the target height in meters assuming a 4/3 effective
    %   %   Earth radius given a range of 300 km, a sensor height of 10 meters,
    %   %   and an elevation angle of 0.5 degrees.
    %   R    = 300e3;
    %   anht = 10;
    %   el   = 0.5;
    %   range2height(R,anht,el)
    %
    %   % Example 2: 
    %   % Calculate the target height assuming a flat Earth, free space
    %   % propagation with a curved Earth, a 4/3 effective Earth radius, and a
    %   % CRPL atmospheric model. Assume a range of 200 km and an antenna 
    %   % height of 100 meters. Plot results over a range of elevation angles 
    %   % from 0 to 5 degrees.
    %   R    = 200e3;
    %   anht = 100;
    %   el   = 0:0.1:5;
    %   
    %   % Calculate target height assuming flat Earth
    %   tgthtFlat = range2height(R,anht,el,'Method','Flat');
    %   
    %   % Calculate target height assuming free space propagation with a curved
    %   % Earth
    %   r0          = physconst('EarthRadius'); 
    %   tgthtFS     = range2height(R,anht,el,'Method','Curved', ...
    %       'EffectiveEarthRadius',r0); 
    % 
    %   % Calculate target height assuming 4/3 effective Earth radius using 
    %   % default values of range2height function. 
    %   tgthtEffRad = range2height(R,anht,el);
    %   
    %   % Calculate target height using CRPL atmospheric model. 
    %   tgthtCRPL   = range2height(R,anht,el,'Method','CRPL');
    %   
    %   % Plot results
    %   figure()
    %   plot(el,[tgthtFlat(:) tgthtFS(:) tgthtEffRad(:)], ...
    %       'LineWidth',1.5)
    %   hold on
    %   plot(el,tgthtCRPL,'--','LineWidth',1.5)
    %   grid on 
    %   xlabel('Elevation Angle (deg)')
    %   ylabel('Target Height (m)')
    %   legend('Flat','Free Space','4/3 Earth','CRPL','Location','Best') 
    %   title('Target Height Estimation')
    %
    %   See also height2range, height2grndrange, el2height, height2el,
    %   refractionexp, radarvcd, blakechart.
    
    %   Copyright 2021-2022 The MathWorks, Inc.  
    
    %   Reference
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off., 1959.
    %   [3] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %       Artech House, 2013.
    '''
    isIterative = False
    if maxNumIter is not None:
        isIterative = True
    
    # Calculate height
    match method.lower():
        case 'curved':
            ht = el2height(el, anht, R, 'Curved', Re)
        case 'flat':
            ht = el2height(el, anht, R, 'Flat')
        case _: # CRPL
            # Convert m to km
            anht = anht * 1e-3 # Antenna height (km)
            R = R * 1e-3 # Target height (km)
            
            # Perform CRPL. Height is relative to antenna. 
            if isIterative:
                ht, _, _ = iterativeCRPL(R, anht, el, Ns, rexp, maxNumIter, tol) # Height (km)
            else:
                ht, _, _ = nonIterativeCRPL(R, anht, el, Ns, rexp) # Height (km)
            
            # Convert km to m
            ht = np.array(ht) * 1e3 # Height (m)
            
    return ht

def iterativeCRPL(Rrow, anhtRow, elRow, Ns, rexp, maxNumIter, tol):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %iterativeCRPL Target height using iterative CRPL method
    %   HT = iterativeCRPL(R,ANHT,EL,NS,REXP,MAXNUMITER,TOL) calculates the
    %   height HT in km using an iterative method assuming a CRPL atmosphere.
    %
    %   The inputs are as follows:
    %      R          = Range (km). M-length row vector.
    %      ANHT       = Antenna height (km). M-length row vector.
    %      EL         = Elevation angle (deg). M-length row vector.
    %      NS         = Refractivity (N-units). Scalar.
    %      REXP       = Decay constant (1/km). Scalar.
    %      MAXNUMITER = Maximum number of iterations. Scalar.
    %      TOL        = Tolerance. Scalar.
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also range2height. 
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''

    # Setup
    Reff = (4 / 3) * Rearth
    
    # Calculate hmax
    hmax = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Rearth) * 1e-3 # Eqn 10
    hmax = hmax - anhtRow # Height relative to antenna
    
    # Calculate hmin
    hmin = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Reff) * 1e-3 # Eqn 10
    hmin = hmin - anhtRow # Height relative to antenna
    
    # hmax ~= hmin ~= 0
    idxNot0 = (hmax > sqrt(eps())) & (hmin > sqrt(eps()))
    if ~np.any(idxNot0):
        # Return if there are no target heights that are not about 0
        h1n = hmin + anhtRow # Height relative to surface 
        numIter = 0
        R1n = Rrow 
        return h1n, numIter, R1n
    
    # Calculate h11
    idxAdjust = abs(hmax - hmin) <= sqrt(eps())
    if np.any(idxAdjust):
        # Adjust hmax/hmin for case where hmax ~ hmin. This typically occurs at
        # angles close to 90 degrees. 
        hmax[idxAdjust] = hmax[idxAdjust] + 0.1
        hmin[idxAdjust] = hmin[idxAdjust] - 0.1
        
    h11 = np.array(hmin)
    h11[idxNot0] = (hmax[idxNot0] + hmin[idxNot0]) / 2 # Eqn 12
    
    # Calculate R11
    elRow = deg2rad(elRow)
    R11 = np.array(Rrow)
    R11[idxNot0] = rangeIntegralCRPL(h11[idxNot0], anhtRow[idxNot0], elRow[idxNot0], Ns, rexp)
    
    # Calculate h12
    h12 = h11
    h12[idxNot0] = h11[idxNot0] + (hmax[idxNot0] - hmin[idxNot0]) / 4 # Eqn 13
    R12 = np.array(Rrow)
    R12[idxNot0] = rangeIntegralCRPL(h12[idxNot0], anhtRow[idxNot0], elRow[idxNot0], Ns, rexp)
    
    # Setup values for first iteration
    R1n = R12
    R1min1 = R12
    R1min2 = R11
    h1min1 = h12
    h1min2 = h11
    h1n = h1min1
    denom = np.array(R1min2 - R1min1)
    
    # Calculate first tolerance
    T = np.array(abs(R12 - Rrow)) # Eqn 15
    idxComplete = np.array(T <= tol)
    
    # Iterate
    numIter = 0
    while (numIter < maxNumIter) and np.any(~idxComplete):
        # Update iteration count
        numIter = numIter + 1
        
        # Calculate new values
        idxIter = ~idxComplete
        denom[idxIter] = R1min2[idxIter] - R1min1[idxIter]
        if abs(denom) < eps(): # Don't let the denominator become zero
            signVal = ones(size(denom))
            signVal[denom < 1] = -1
            denom = eps() * signVal

        h1n[idxIter] = (Rrow[idxIter] - R1min1[idxIter]) * (h1min2[idxIter] - h1min1[idxIter]) / (denom[idxIter]) + h1min1[idxIter] # Eqn 14
        R1n[idxIter] = rangeIntegralCRPL(h1n[idxIter], anhtRow[idxIter], elRow[idxIter], Ns, rexp) # Eqn 16
        
        # Calculate threshold
        T[idxIter] = abs(R1n[idxIter] - Rrow[idxIter]) # Eqn 15
        
        # Update values for next iteration
        R1min2[idxIter] = R1min1[idxIter]
        R1min1[idxIter] = R1n[idxIter]
        h1min2[idxIter] = h1min1[idxIter]
        h1min1[idxIter] = h1n[idxIter]
        
        # Set output for values that satisfy tolerance
        idxComplete[idxIter] = T[idxIter] <= tol 

    h1n = h1n + anhtRow # Height relative to surface

    return h1n, numIter, R1n

def nonIterativeCRPL(Rrow, anhtRow, elRow, Ns, rexp):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %nonIterativeCRPL Target height using non-iterative CRPL method
    %   HT = nonIterativeCRPL(R,ANHT,EL,NS,REXP) calculates the height HT in km
    %   using a non-iterative method assuming a CRPL atmosphere.
    %
    %   The inputs are as follows:
    %      R          = Range (km). M-length row vector.
    %      ANHT       = Antenna height (km). M-length row vector.
    %      EL         = Elevation angle (deg). M-length row vector.
    %      NS         = Refractivity (N-units). Scalar. 
    %      REXP       = Decay constant (1/km). Scalar. 
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also range2height. 
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''    
    # Setup 
    Reff = (4 / 3) * Rearth
    
    # Calculate hmin
    hmin = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Reff) * 1e-3 # Eqn 10
    hmin = hmin - anhtRow # Height relative to antenna
    
    # Calculate Rhmin
    elRow = deg2rad(elRow)
    Rhmin = rangeIntegralCRPL(hmin, anhtRow, elRow, Ns, rexp) # Eqn 2 using hmin
    
    # Calculate delR
    delR = abs(Rrow - Rhmin) # Eqn 43
    
    # Calculate r1, theta1, and rho1
    r0 = Rearth * 1e-3
    r1 = r0 + hmin # Eqn 39
    rho0 = Ns * 1e-6
    theta1 = acos((1 + rho0) * cos(elRow) / ((1 + rho0 * exp(-rexp * hmin)) * (1 + hmin / r0))) # Eqn 40
    rho1 = rho0 * exp(-rexp * hmin) # Eqn 41
    
    # Calculate b1 and c1
    gammaConst = rho1 * rexp / (1 + rho1) # Eqn 23
    b1 = 2 * ((2 * gammaConst - 1 / r1) * sin(theta1) ** 2 - gammaConst + 1 / r1) # Eqn 26
    c1 = sin(theta1) ** 2 # Eqn 27
    
    # Calculate h1 
    h1 = b1 * (delR / (2 * (1 + rho1))) ** 2 + delR * sqrt(c1) / (1 + rho1) + hmin # Eqn 44
    h1 = h1 + anhtRow # Height relative to surface
    
    return h1

def el2height(
    el,
    anht,
    SR,
    model = 'Curved',
    Re = (4 / 3) * Rearth,
):
    '''
    %el2height Convert target initial elevation angle to height
    %   HT = el2height(EL,ANHT,SR) returns the target height in meters. EL is
    %   the target elevation angle in degrees, ANHT is the sensor height in
    %   meters, and SR is the slant range between the target and the sensor in
    %   meters. The computation assumes a curved Earth model with 4/3 effective
    %   Earth radius, which is an approximation used for modeling refraction
    %   effects in the troposphere. Heights are referenced to the ground.
    %
    %   EL, ANHT, and SR can either be scalars or M-length vectors.
    %
    %   HT = el2height(EL,ANHT,SR,MODEL) specifies the Earth model used to
    %   compute the target height as one of 'Flat' | 'Curved', where the
    %   default is 'Curved'.
    %
    %   HT = el2height(EL,ANHT,SR,MODEL,RE) specifies the effective Earth
    %   radius in meters as a positive scalar, RE, where the default is 4/3 of
    %   the Earth radius. The effective Earth radius is an approximation used
    %   for modeling refraction effects in the troposphere. This input is
    %   ignored when you set Model to 'Flat'.
    %
    %   % Example:
    %   %   Determine the target height in meters given an elevation angle of 
    %   %   0.5 degrees, a sensor height of 10 meters, and a slant range of 300
    %   %   km. 
    %   el   = 0.5;
    %   anht = 10;
    %   SR   = 300e3;
    %   el2height(el,anht,SR)
    %
    %   See also height2el, horizonrange, depressionang, grazingang,
    %   effearthradius.
    
    %   Copyright 2020-2022 The MathWorks, Inc.  
    
    %   Reference
    %   [1] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %       Artech House, 2013.
    '''
    
    if model[0] == 'C':
        # Curved Earth calculation
        val = (Re + anht) ** 2 + SR ** 2 + 2 * (Re + anht) * SR * sind(el) 
        ht = sqrt(val) - Re # Equation 7.9 in Barton (m)
    
        # Assume flat Earth calculation for when el ~ 90
        idx90 = 90 - abs(el) < 1
        htFlat = anht + SR * sind(el) # (m)

        if numel(ht) > 1:
            ht[idx90] = htFlat[idx90]
        elif idx90:
            ht = htFlat
    else:
        # Flat Earth calculation
        ht = anht + SR * sind(el) # (m)
        idxNaN = isnan(ht)
        if numel(anht) > 1:
            ht[idxNaN] = anht[idxNaN] # Edge case when R*sind(theta) = inf*0 = NaN
        else:
            ht[idxNaN] = anht # Edge case when R*sind(theta) = inf*0 = NaN
    
    return ht

def cart2pol(x, y, z = None):
    '''
    %CART2POL Transform Cartesian to polar coordinates.
    %   [TH,R] = CART2POL(X,Y) transforms corresponding elements of data stored
    %   in Cartesian coordinates X,Y to polar coordinates (angle TH and radius
    %   R). The arrays X and Y must have compatible sizes. In the simplest
    %   cases, they can be the same size or one can be a scalar. Two inputs
    %   have compatible sizes if, for every dimension, the dimension sizes of
    %   the inputs are either the same or one of them is 1. TH is returned in
    %   radians.
    %
    %   [TH,R,Z] = CART2POL(X,Y,Z) transforms corresponding elements of data
    %   stored in Cartesian coordinates X,Y,Z to cylindrical coordinates (angle
    %   TH, radius R, and height Z). The arrays X,Y, and Z must have compatible
    %   sizes. TH is returned in radians.
    %
    %   Class support for inputs X,Y,Z:
    %      float: double, single
    %
    %   See also CART2SPH, SPH2CART, POL2CART.
    
    %   Copyright 1984-2021 The MathWorks, Inc. 
    '''
    th = atan2(y, x)
    r = hypot(x, y)

    if z is not None:
        return th, r, z
    else:
        return th, r

def pol2cart(th, r, z = None):
    '''
    %POL2CART Transform polar to Cartesian coordinates.
    %   [X,Y] = POL2CART(TH,R) transforms corresponding elements of data stored
    %   in polar coordinates (angle TH, radius R) to Cartesian coordinates X,Y.
    %   The arrays TH and R must have compatible sizes. In the simplest cases,
    %   they can be the same size or one can be a scalar. Two inputs have
    %   compatible sizes if, for every dimension, the dimension sizes of the
    %   inputs are either the same or one of them is 1. TH must be in radians.
    %
    %   [X,Y,Z] = POL2CART(TH,R,Z) transforms corresponding elements of data
    %   stored in cylindrical coordinates (angle TH, radius R, height Z) to
    %   Cartesian coordinates X,Y,Z. The arrays TH, R, and Z must have
    %   compatible sizes.  TH must be in radians.
    %
    %   Class support for inputs TH,R,Z:
    %      float: double, single
    %
    %   See also CART2SPH, CART2POL, SPH2CART.
    
    %   L. Shure, 4-20-92.
    %   Copyright 1984-2021 The MathWorks, Inc. 
    '''
    x = r * cos(th)
    y = r * sin(th)

    if z is not None:
        return x, y, z
    else:
        return x, y
    
# Plotting helper functions
# -----------------------------------------+
def plotCustomAxes(hAxes, rngX, rngY, htX, htY, angX, angY, curveColor):
    # Plot constant height curves
    hHeightAxes = line(hAxes, htX[:-1, :], htY[:-1, :], HandleVis = 'off', HitTest = 'off', Tag = 'ConstantHeightCurves')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hHeightAxes, 'Color', curveColor)
    # numH = numel(hHeightAxes)
    # arrayfun(@(idx) set(hHeightAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    # Plot constant range curves
    hRangeAxes = line(hAxes,rngX, rngY, HandleVis = 'off', HitTest = 'off', Tag = 'ConstantRangeCurves')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hRangeAxes, 'Color', curveColor)
    # numH = numel(hRangeAxes)
    # arrayfun(@(idx) set(hRangeAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    # Plot angle pokes
    hAngleAxes = line(hAxes, angX, angY, HandleVis = 'off', HitTest = 'off', Tag = 'ConstantAngleLines')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngleAxes, 'Color', curveColor)
    # numH = numel(hAngleAxes)
    # arrayfun(@(idx) set(hAngleAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    hAngleAxesEnd = line(hAxes, angX[:, -1], angY[:, -1], HandleVis = 'off', HitTest = 'off', Tag = 'ConstantAngleEndPoints')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngleAxesEnd, 'Color', '--mw-graphics-borderColor-axes-primary')
    # numH = numel(hAngleAxesEnd)
    # arrayfun(@(idx) set(hAngleAxesEnd(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
# function fixOverlappingText(hText)
#     % Fix overlapping text for angular axes
#     numH = numel(hText);
#     idxDelete = false(1,numH);
#     iiprev = 1;
#     if numH > 1
#         for ii = 2:numH
#             x2 = hText(ii).Extent(1);
#             y2 = hText(ii).Extent(2);
#             w2 = hText(ii).Extent(3);
#             h2 = hText(ii).Extent(4);
#             rectPts = [x2 y2; x2 (y2 + h2); (x2 + w2) (y2 + h2); x2 (y2 + h2)];
            
#             % Check for overlapping x
#             minX = hText(iiprev).Extent(1);
#             maxX = hText(iiprev).Extent(1) + hText(iiprev).Extent(3);
#             idxOverlapX = any(rectPts(:,1) >= minX & rectPts(:,1) <= maxX);
            
#             % Check for overlapping y
#             minY = hText(iiprev).Extent(2);
#             maxY = hText(iiprev).Extent(2) + hText(iiprev).Extent(4);
#             idxOverlapY = any(rectPts(:,2) >= minY & rectPts(:,2) <= maxY);
            
#             delCond = idxOverlapX && idxOverlapY;
#             [idxDelete,iiprev] = updateDeleteIdx(delCond,idxDelete,iiprev,ii);
#         end
#         delete(hText(idxDelete));
#     end
#     end
    
# function [ticks,tickLabel] = fixOverlappingTicks(hAxes,hText,ticks,tickLabel,type)
#     % Estimate text dimensions
#     switch type
#         case 'x'
#             % Estimate tick text width
#             tickLength = hText(1).Extent(3);
#             numChar = numel(hText(1).String);
#             tickLength = tickLength/numChar;
            
#             % Estimate tick extent in horizontal direction
#             formatType = hAxes.XAxis.TickLabelFormat;
#             xTickLabelsChar = arrayfun(@(x) num2str(x,formatType),tickLabel,'UniformOutput',false);
#             numChars = cellfun(@(x) numel(x) + 1,xTickLabelsChar); % Add extra character for spacing
#             halfTickLength = (tickLength.*numChars(:))./2;
#         otherwise
#             % Estimate tick text height
#             tickLength = hText(1).Extent(4);
#             halfTickLength = (tickLength)/2;
#     end
    
#     % Calculate text extent
#     textExtent = [(ticks(:) - halfTickLength) (ticks(:) + halfTickLength)]; % Min/max y2
    
#     % Identify overlapping tick marks
#     numTicks = numel(ticks);
#     idxDelete = false(1,numTicks);
#     iiprev = 1;
#     if numTicks > 1
#         for ii = 2:numTicks
#             delCond = any((textExtent(ii,:) >= textExtent(iiprev,1)) ...
#                 & (textExtent(ii,:) <= textExtent(iiprev,2)));
#             [idxDelete,iiprev] = updateDeleteIdx(delCond,idxDelete,iiprev,ii);
#         end
#     end
    
#     % Delete overlapping tick marks
#     ticks(idxDelete) = [];
#     tickLabel(idxDelete) = [];
#     end
    
# function [idxDelete,iiprev] = updateDeleteIdx(delCond,idxDelete,iiprev,ii)
#     % Update indices to delete
#     if delCond
#         % Add index to list of those to be deleted. Do not update previous
#         % index.
#         idxDelete(ii) = true;
#     else
#         % Update previous index
#         iiprev = ii;
#     end
#     end
    
# function addCustomDataTips(hPlot,labels,vals)
#     % Add custom data tips (more than 2 data tip rows are allowed) to the plots
#     % in hPlot.
#     numPlots = numel(hPlot);
#     numDTs = numel(labels);
#     for ii = 1:numPlots
#         for jj = 1:numDTs
#             if jj <= 2
#                 hPlot(ii).DataTipTemplate.DataTipRows(jj).Label = labels{jj};
#                 hPlot(ii).DataTipTemplate.DataTipRows(jj).Value = vals{jj};
#             else
#                 % For additional data tip rows past 2
#                 row = dataTipTextRow(labels{jj},vals{jj});
#                 hPlot(ii).DataTipTemplate.DataTipRows(end+1) = row;
#             end
#         end
#     end
#     end
    
# function     idx = setDataTipIndex(vcpX,vcpY,xlimMax,ylimMax)
#     % Try to find a datatip index within the x/y limits of the figure
#     idx = find((vcpX < xlimMax) & (vcpY < ylimMax),1,'first');
    
#     % Otherwise, just set it to the first index
#     if isempty(idx)
#         idx = 1;
#     end
#     end
    
# function applyPatchColor(hPatch,fc,idx)
#     if iscell(fc)
#         % Custom colors 
#         hPatch.FaceColor = fc{idx}; 
#     else
#         % Default colors
#         thisColorIndex = fc(idx); 
#         patchColor = sprintf('--mw-graphics-colorOrder-%d-primary',thisColorIndex); 
#         matlab.graphics.internal.themes.specifyThemePropertyMappings(...
#             hPatch,'FaceColor',patchColor);
    
#     end
#     end
    
# function applyPlotColor(hPlot,ec,idx)
#     if iscell(ec)
#         % Custom colors
#         hPlot.Color = ec{idx}; 
#     else
#         % Default colors
#         thisColorIndex = ec(idx);
#         lineColor = sprintf('--mw-graphics-colorOrder-%d-primary',thisColorIndex); 
#         matlab.graphics.internal.themes.specifyThemePropertyMappings(...
#             hPlot,'Color',lineColor);
#     end
#     end

if __name__ == "__main__":
    pass