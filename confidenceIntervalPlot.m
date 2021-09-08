classdef confidenceIntervalPlot < matlab.graphics.chartcontainer.ChartContainer & ...
        matlab.graphics.chartcontainer.mixin.Legend
    %confidenceIntervalPlot Create a mean line within a shaded confidence interval area.
    %   confidenceIntervalPlot(x,y) create a line which passes through the 
    %   means of the y-values for each unique x-value. Plot this line within
    %   a shaded area covering a 95% confidence interval for each 
    %   unique x-value. x and y must be numeric vectors of equal length.
    %
    %   confidenceIntervalPlot(x,y,alpha) create a line which passes through 
    %   the means of the y-values for each unique x-value. Plot this line 
    %   within a shaded area covering a 100 * (1 - alpha)% confidence interval 
    %   for each unique x-value. x and y must be numeric vectors of equal 
    %   length.
    %
    %   confidenceIntervalPlot() create an empty confidence interval plot.
    %
    %   confidenceIntervalPlot(___,Name,Value) specifies additional options
    %   for the confidence interval plot using one or more name-value pair
    %   arguments. Specify the options after all other input arguments.
    %
    %   confidenceIntervalPlot(parent,___) creates the confidence interval 
    %   plot in the specified parent.
    %
    %   h = confidenceIntervalPlot(___) returns the confidenceIntervalPlot
    %   object. Use h to modify properties of the plot after creating it.
    
    %   Copyright 2021 The MathWorks, Inc.

    properties
        XData (1,:) {mustBeNumeric} = []
        YData (1,:) {mustBeNumeric} = []

        % y values used to plot the center line
        CenterYData (1,:) {mustBeNumeric} = []

        % Confidence level for the confidence interval about the mean line
        Alpha (1,1) {mustBeScalarOrEmpty} = 0.05;

        % Upper and lower bound data of the shaded area encapsulating the 
        % mean line
        UpperBoundData (1,:) {mustBeNumeric} = []
        LowerBoundData (1,:) {mustBeNumeric} = []

        % Bins to group together x-values rather than use unique x-values
        Edges (1,:) {mustBeNumeric} = []
 
        % Title and subtitle of the confidence interval plot
        TitleText (1,:) char = ''
        SubtitleText (:,:) char = ''

        % Color and transparency of the shaded area surrounding the center
        % line
        ShadeColor (1,3) {mustBeNumeric} = [0 0 1]
        ShadeAlpha (1,1) {mustBeScalarOrEmpty} = 0.1

        % Color and width of the center line
        CenterLineColor (1,3) {mustBeNumeric} = [0 0 1]
        CenterLineWidth (1,1) {mustBeScalarOrEmpty} = 0.5

        % Color and width of the border lines around the shaded area about
        % the center line
        BorderLinesColor (1,3) {mustBeNumeric} = [1 0 0]
        BorderLinesWidth (1,1) {mustBeScalarOrEmpty} = 0.5

        % Properties of the raw data markers if ShowRawData is true
        RawDataMarker char {mustBeMarker} = 'o'
        RawDataMarkerColor (1,3) {mustBeNumeric} = [0 0 1]
        RawDataMarkerSize (1,1) {mustBeScalarOrEmpty} = 10;

        % Properties of the center data markers if ShowCenterData is true
        CenterDataMarker char {mustBeMarker} = '*'
        CenterDataMarkerColor (1,3) {mustBeNumeric} = [0 0 1]
        CenterDataMarkerSize (1,1) {mustBeScalarOrEmpty} = 10

        % Number of points between each pair of points in CenterXData, 
        % determines smoothness of curves via interpolation
        NumSteps (1,1) {mustBeScalarOrEmpty, mustBeNonnegative} = 0

        % Boolean indicating whether to show the raw data with scatter
        ShowRawData (1,1) matlab.lang.OnOffSwitchState {mustBeScalarOrEmpty} = true

        % Boolean indicating whether to show the data used to create the
        % center line
        ShowCenterData (1,1) matlab.lang.OnOffSwitchState {mustBeScalarOrEmpty} = false

        % String describing the method by which the bounds of the shaded
        % area are determined (confidence interval or manual)
        BoundDataMode (1,:) char {mustBeAutoManual} = 'auto'

        % String describing the method by which the bounds of the shaded
        % area are determined (confidence interval or manual)
        CenterYDataMode (1,:) char {mustBeAutoManual} = 'auto'
    end

    properties (SetAccess = private)
        % read-only property for x values used to plot the center line
        CenterXData (1,:) {mustBeNumeric} = []
    end

    properties (Access = protected)
        % Used for saving to .fig files
        ChartState = []
    end

    properties(Access = private,Transient,NonCopyable)
        CenterLine (1,:) matlab.graphics.primitive.Line
        LowerBoundLine (1,:) matlab.graphics.primitive.Line
        UpperBoundLine (1,:) matlab.graphics.primitive.Line
        ShadedArea (1,:) matlab.graphics.primitive.Patch
        RawDataScatter (1,:) matlab.graphics.chart.primitive.Scatter
        CenterDataScatter (1,:) matlab.graphics.chart.primitive.Scatter
    end

    methods
        function obj = confidenceIntervalPlot(varargin)
            % Initialize list of arguments
            args = varargin;
            leadingArgs = cell(0);

            % Check if the first input argument is a graphics object to use as parent.
            if ~isempty(args) && isa(args{1},'matlab.graphics.Graphics')
                % confidenceIntervalPlot(parent, ___)
                leadingArgs = args(1);
                args = args(2:end);
            end

            % Check for optional positional arguments.
            if ~isempty(args) && isnumeric(args{1})
                if numel(args) >= 2 && mod(numel(args), 2) == 0
                    % confidenceIntervalPlot(x,y)
                    % confidenceIntervalPlot(x,y,Name,Value)
                    x = args{1};
                    y = args{2};

                    leadingArgs = [leadingArgs {'XData', x, 'YData', y}];
                    args = args(3:end);
                elseif numel(args) >= 3 && mod(numel(args), 2) == 1 && isnumeric(args{3})
                    % confidenceIntervalPlot(x,y,alpha)
                    % confidenceIntervalPlot(x,y,alpha,Name,Value)
                    x = args{1};
                    y = args{2};
                    alpha = args{3};

                    leadingArgs = [leadingArgs {'XData', x, 'YData', y, 'Alpha', alpha}];
                    args = args(4:end);
                else
                    error('Invalid number of input arguments.');
                end

            end

            % Combine positional arguments with name/value pairs.
            args = [leadingArgs args];

            % Call superclass constructor method
            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
        end

    end

    methods(Access = protected)

        function setup(obj)
            % Create the axes
            ax = getAxes(obj);

            % Remove axes toolbar
            ax.Toolbar = [];
            
            % Call the load method in case of loading from a fig file
            loadstate(obj);

            % Create graphics objects
            obj.ShadedArea = patch(ax, 'LineStyle', 'None', 'DisplayName', 'Bounded Area');
            obj.CenterLine = line(ax, 1, 1, 'DisplayName', 'Center Line');
            obj.LowerBoundLine = line(ax, 1, 1, 'DisplayName', 'Lower Bound Line');
            obj.UpperBoundLine = line(ax, 1, 1, 'DisplayName', 'Upper Bound Line');

            ax.NextPlot = 'add';
            obj.RawDataScatter = scatter(ax, obj.XData, obj.YData, ...
                'DisplayName', 'Raw Data');
            obj.CenterDataScatter = scatter(ax, obj.CenterXData, obj.CenterYData, ...
                'DisplayName', 'Center Line Data');
        end

        function update(obj)
            ax = getAxes(obj);

            % Update title
            title(ax, obj.TitleText, obj.SubtitleText);
            
            % Throw an error if the x and y data are not the same size
            if numel(obj.XData) ~= numel(obj.YData)
                error('X-data and Y-data must be the same size.');
            end

            % If there is no XData, there is no need to fill the axes
            if isempty(obj.XData)
                obj.ShadedArea.Visible = 'off';
                return
            else
                % Make patch visible if there is nonempty XData (passing return
                % above)
                obj.ShadedArea.Visible = 'on';
            end

            % Preprocessing: store valid x and y data which exclude all x 
            % and y values in the raw data which have NaN x-value
            validXIndices = ~isnan(obj.XData);
            validYData = obj.YData(validXIndices);
            validXData = obj.XData(validXIndices);

            % After this if statement, iCenterXData is an index vector,
            % the same length as the raw data (in-bin if applicable), that 
            % indicates which bin/unique x-value group each data point belongs to
            if isempty(obj.Edges)
                % Find unique x data and index of each element of XData in
                % UniqueXData
                [obj.CenterXData, ~, CenterXDataIndices] = unique(obj.XData);
            else
                % Group y-values according to bins in edges
                CenterXDataIndices = discretize(obj.XData, obj.Edges);

                % Get indices which are not NaN (i.e. x data indices for data
                % which fall in the specified bins)
                validBinIndices = ~isnan(CenterXDataIndices);

                % Handling NaN values: for bin indices which are NaN, discard
                % corresponding x and y values
                CenterXDataIndices = CenterXDataIndices(validBinIndices);
                validYData = validYData(validBinIndices);
                validXData = validXData(validBinIndices);

                % Compute the center x-values to be the mean of the values 
                % in each bin
                obj.CenterXData = accumarray(CenterXDataIndices(:), validXData(:), [], @mean);
            end

            % Number of elements in CenterXData
            numCenterXData = numel(obj.CenterXData);

            % Contains a 1 if covered by a patch (at least 1 non-NaN y-value), 
            % otherwise 0 (gap due to all NaN y-values).
            patchIndices = false(1, numCenterXData);
            
            % Maximum number of vertices on a face
            maxFaceVertices = 0;

            % Total number of faces
            numFaces = 0;

            % Whether the last center x-value is covered by a patch
            lastPatchIdx = 0;

            % Index of the last center x-value which is at the start of a patch
            lastPatchStartIdx = 0;

            % Set iCenterXData to NaN where y-values are NaN. Then below,
            % elements of iCenterXData will never equal centerXIdx for
            % invalid YData.
            CenterXDataIndices(isnan(obj.YData)) = NaN;

            % logical vector for indices which have all NaNs
            yAllNaNIndices = false(1, numCenterXData);

            for centerXIdx = 1:numCenterXData
                % Get all y-values associated with the current center
                % x-value
                yVals = validYData(CenterXDataIndices == centerXIdx);

                % If there is at least 1 non-NaN y-value, save the
                % index as an index covered by a patch. Otherwise save the 
                % index as uncovered by a patch.
                if ~isempty(yVals)
                    patchIndices(centerXIdx) = 1;

                    % If we just reached a new patch, save index as the start
                    % of the last patch 
                    if lastPatchIdx == 0
                        lastPatchStartIdx = centerXIdx;
                    end
                else
                    patchIndices(centerXIdx) = 0;
                    yAllNaNIndices(centerXIdx) = true;

                    % If we've just finished a patch, add 1 to the total
                    % number of faces and update the maximum number of
                    % vertices.
                    if lastPatchIdx == 1
                        numFaces = numFaces + 1;
                        lastPatchLength = centerXIdx - lastPatchStartIdx + 1;
                        maxFaceVertices = max([maxFaceVertices, 2 * lastPatchLength]);
                    end
                end

                lastPatchIdx = patchIndices(centerXIdx);

            end

            % If the last center x-value is covered by a patch, add 1 to
            % the total number of faces and update the maximum number of
            % vertices.
            if ~isempty(patchIndices) && patchIndices(end) == 1
                numFaces = numFaces + 1;
                lastPatchLength = numel(patchIndices) - lastPatchStartIdx + 1;
                maxFaceVertices = max([maxFaceVertices, 2 * lastPatchLength]);
            end

            % Update validXData, iCenterXData, and validYData so that they 
            % exclude points with nan y-value
            validXData = validXData(~isnan(validYData));
            CenterXDataIndices = CenterXDataIndices(~isnan(validYData));
            validYData = validYData(~isnan(validYData));

            % Check the center data mode to determine whether we need to
            % compute the center y-data
            if strcmp(obj.CenterYDataMode, 'auto')
                % Remove the NaN values from iCenterXData before
                % calculating
                CenterXDataIndices = CenterXDataIndices(~isnan(CenterXDataIndices));
                centerYDataCol = accumarray(CenterXDataIndices(:), validYData(:), [], @mean);
                obj.CenterYData = centerYDataCol(:)';
            end

            % For x-values with all nan y-values, set center y-value to
            % NaN--then plotting lines will result in the correct behavior
            obj.CenterYData(yAllNaNIndices) = NaN;

            % Check the bound data mode to determine whether we need to
            % compute the upper and lower bound data
            if strcmp(obj.BoundDataMode, 'auto')
                % For each group of indices corresponding to a unique x value,
                % find the mean of the corresponding y values to get the mean
                % data points. Repeat for standard deviation
                StdDevData = accumarray(CenterXDataIndices(:), validYData(:), [], @std);
                StdDevData = StdDevData.';
    
                % To get the appropriate confidence level, convert the 
                % confidence level to z-score
                z_score = norminv(1 - obj.Alpha);
    
                % Get lower bound and upper bound data
                obj.LowerBoundData = obj.CenterYData - StdDevData * z_score;
                obj.UpperBoundData = obj.CenterYData + StdDevData * z_score;
            end

            % By default, assume no interpolation
            finerCenterXData = obj.CenterXData;
            finerCenterYData = obj.CenterYData;
            finerLowerBoundData = obj.LowerBoundData;
            finerUpperBoundData = obj.UpperBoundData;

            % Interpolate all lines/patches
            if obj.NumSteps ~= 0

                % total number of points (with interpolation) used to plot the center line
                totalNumPoints = (numCenterXData - 1) * obj.NumSteps + numCenterXData;

                finerCenterXData = zeros(1, totalNumPoints);

                % For each gap between the original data points, add
                % additional points such that there are NumSteps points
                % between each pair of the original data points.
                for i = 1:(numCenterXData - 1)
                    startIdx = 1 + (i - 1) * (obj.NumSteps + 1);
                    endIdx = 1 + i * (obj.NumSteps + 1);
                    finerCenterXData(startIdx:endIdx) = linspace(obj.CenterXData(i), ...
                        obj.CenterXData(i + 1), obj.NumSteps + 2);
                end
                
                % Interpolate the two bounds, the center line, and the
                % shaded area/patch data, and update the coordinates
                finerLowerBoundData = interp1(obj.CenterXData, obj.LowerBoundData, ...
                    finerCenterXData, 'spline');
                finerUpperBoundData = interp1(obj.CenterXData, obj.UpperBoundData, ...
                    finerCenterXData, 'spline');
                finerCenterYData = interp1(obj.CenterXData, obj.CenterYData, ...
                    finerCenterXData, 'spline');
            end

            % Boolean indicating whether the last index belonged to a patch
            started_patch = false;

            % Faces and vertices for drawing final patch
            finerMaxFaceVertices = maxFaceVertices + 2 * obj.NumSteps * (maxFaceVertices / 2 - 1);
            faces = NaN(numFaces, finerMaxFaceVertices);
            vertices = [];

            % Index of the last face added to faces and the last vertex
            % added to vertices
            vertexIdx = 1;
            faceIdx = 1;

            for centerXIdx = 1:numCenterXData
                
                if patchIndices(centerXIdx) == 0 || centerXIdx == numCenterXData
                    % If the current x-value shouldn't be covered by a patch
                    % and the last index belonged to a patch, update faces 
                    % and vertices to include the corresponding patch

                    if started_patch
                        % If the last x-value should be covered by a patch,
                        % set the end index of the patch to be the current
                        % index.
                        if centerXIdx == numCenterXData
                            endIdx = centerXIdx;
                        else
                            endIdx = centerXIdx - 1;
                        end

                        % Convert the indices of the x-values in
                        % CenterXData to the indices of the same x-values
                        % in finerCenterXData
                        startIdx = (startIdx - 1) * (obj.NumSteps + 1) + 1;
                        endIdx = (endIdx - 1) * (obj.NumSteps + 1) + 1;
                        centerIndices = startIdx:endIdx;

                        % Identify vertices of the polygon to draw for the shaded area
                        % and draw it
                        patchXData = [finerCenterXData(centerIndices), fliplr(finerCenterXData(centerIndices))];
                        patchYData = [finerLowerBoundData(centerIndices), fliplr(finerUpperBoundData(centerIndices))];

                        % Add a new row to faces and add all vertices for
                        % this face to vertices
                        endVertexIdx = vertexIdx + length(patchXData) - 1;
                        newFace = vertexIdx:endVertexIdx;
                        faces(faceIdx, 1:numel(newFace)) = newFace;
                        vertices = [vertices; patchXData(:) patchYData(:)]; %#ok<AGROW> 

                        % Update indices for faces and vertices
                        faceIdx = faceIdx + 1;
                        vertexIdx = endVertexIdx + 1;
                        started_patch = false;
                    end
                else
                    % If the current x-value should be covered by a patch
                    % and the last index didn't belong to a patch, update
                    % the flag and save the x-value index
                    if ~started_patch
                        startIdx = centerXIdx;
                        started_patch = true;
                    end
                end
            end

            if numFaces
                obj.ShadedArea.Faces = faces;
                obj.ShadedArea.Vertices = vertices;
                obj.ShadedArea.FaceColor = obj.ShadeColor;
                obj.ShadedArea.FaceAlpha = obj.ShadeAlpha;
            end

            % Update the center and upper/lower bound lines
            Names = {'XData', 'YData', 'Color', 'LineWidth'};
            CenterLineValues = {finerCenterXData, finerCenterYData, ...
                obj.CenterLineColor, obj.CenterLineWidth};
            LowerBoundLineValues = {finerCenterXData, finerLowerBoundData, ...
                obj.BorderLinesColor, obj.BorderLinesWidth};
            UpperBoundLineValues = {finerCenterXData, finerUpperBoundData, ...
                obj.BorderLinesColor, obj.BorderLinesWidth};

            set(obj.CenterLine, Names, CenterLineValues);
            set(obj.LowerBoundLine, Names, LowerBoundLineValues);
            set(obj.UpperBoundLine, Names, UpperBoundLineValues);

            % Plot the data as a scatter plot
            if obj.ShowRawData
                obj.RawDataScatter.XData = validXData;
                obj.RawDataScatter.YData = validYData;
                obj.RawDataScatter.Marker = obj.RawDataMarker;
                obj.RawDataScatter.SizeData = obj.RawDataMarkerSize;
                obj.RawDataScatter.MarkerEdgeColor = obj.RawDataMarkerColor;
                obj.RawDataScatter.Visible = 'on';
            else
                obj.RawDataScatter.Visible = 'off';
            end

            % Plot the data used to create the center line
            if obj.ShowCenterData
                obj.CenterDataScatter.XData = obj.CenterXData;
                obj.CenterDataScatter.YData = obj.CenterYData;
                obj.CenterDataScatter.Marker = obj.CenterDataMarker;
                obj.CenterDataScatter.SizeData = obj.CenterDataMarkerSize;
                obj.CenterDataScatter.MarkerEdgeColor = obj.CenterDataMarkerColor;
                obj.CenterDataScatter.Visible = 'on';
            else
                obj.CenterDataScatter.Visible = 'off';
            end
        end

    end

    methods
        % Title method. Called in update with the TitleText property so
        % that the user can specify the title and subtitle of the Chart.
        function title(obj,txt,subtxt)
            if nargin>=2
                if isnumeric(txt)
                    txt=num2str(txt);
                end
                obj.TitleText = txt;
            end

            if nargin == 3
                if isnumeric(subtxt)
                    subtxt=num2str(subtxt);
                end
                obj.SubtitleText = subtxt;
            end
        end

        function subtitle(obj,subtxt)
            if isnumeric(subtxt)
                subtxt=num2str(subtxt);
            end
            obj.SubtitleText = subtxt;
        end

        function data = get.ChartState(obj)
            % This method gets called when a .fig file is saved
            isLoadedStateAvailable = ~isempty(obj.ChartState);

            if isLoadedStateAvailable
                data = obj.ChartState;
            else
                data = struct;
                ax = getAxes(obj);

                % Get axis limits only if mode is manual.
                if strcmp(ax.XLimMode,'manual')
                    data.XLim = ax.XLim;
                end
                if strcmp(ax.YLimMode,'manual')
                    data.YLim = ax.YLim;
                end
            end
        end

        function loadstate(obj)
            % Call this method from setup to handle loading of .fig files
            data=obj.ChartState;
            ax = getAxes(obj);

            % Look for states that changed
            if isfield(data, 'XLim')
                ax.XLim=data.XLim;
            end
            if isfield(data, 'YLim')
                ax.YLim=data.YLim;
            end
        end
    end
end


function mustBeAutoManual(mode)
    mustBeMember(mode, {'auto','manual'})
end

function mustBeMarker(marker)
    mustBeMember(marker, {'+' , 'o' , ...
                    '*' , '.' , 'x' , 'square' , 'diamond' , ...
                    'v' , '^' , '>' , '<' , 'pentagram' , ...
                    'hexagram' , '|' , '_' , 'none'})
end