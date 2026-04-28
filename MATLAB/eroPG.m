classdef eroPG < handle
    % eroPG v8.0 - Command History (Up/Down Arrow Support)
    %
    % Usage:
    %   pg = eroPG(A);
    %
    % New Feature:
    %   - Click the Command Input box.
    %   - Press UP ARROW to recall previous commands.
    %   - Press DOWN ARROW to navigate back.
    
    properties (SetAccess = private)
        % Backend Data
        OriginalMatrix
        CurrentMatrix
        LastE           
        TotalE          
        
        History = {}
        UndoStack = {}
        
        % --- NEW: Command History Properties ---
        CmdHistory = string.empty; % Stores valid commands executed
        CmdHistoryIdx = 1;         % Points to the current position in history
        
        % GUI Components
        Fig
        MatrixDisplayA      
        MatrixDisplayE      
        MatrixDisplayTotal  
        HistoryList
        
        % Inputs
        CmdInput
        SubVarInput
        SubValInput
    end
    
    methods
        % --- CONSTRUCTOR ---
        function obj = eroPG(A)
            if nargin < 1
                A = sym([1 2 3; 4 5 6; 7 8 9]); 
            end
            
            obj.OriginalMatrix = sym(A);
            obj.CurrentMatrix = obj.OriginalMatrix;
            
            n = size(A, 1);
            obj.LastE = eye(n, 'sym');
            obj.TotalE = eye(n, 'sym');
            
            obj.buildGUI();
            obj.updateDisplay();
        end
    end
    
    methods (Access = private)
        % --- GUI LAYOUT ---
        function buildGUI(obj)
            % 1. Main Figure
            % Added WindowKeyPressFcn to listen for Up/Down keys
            obj.Fig = uifigure('Name', 'ERO Playground v8.0', ...
                'Position', [50, 50, 1150, 650], ...
                'WindowKeyPressFcn', @(src, event) obj.onWindowKeyPress(event));
            
            % 2. Matrix Display Area
            container = uipanel(obj.Fig, ...
                'Position', [20, 250, 850, 380], 'BorderType', 'none');
            
            gridDisp = uigridlayout(container, [1, 2]);
            gridDisp.ColumnWidth = {'1.8x', '1x'}; 
            gridDisp.Padding = [0 0 0 0];
            
            % Panel A (Left)
            panelA = uipanel(gridDisp, 'Title', 'Current Matrix (A)', 'FontSize', 14);
            axA = uiaxes(panelA, 'Position', [10, 10, 520, 330]);
            axis(axA, 'off');
            obj.MatrixDisplayA = text(axA, 0.5, 0.5, '', 'Interpreter', 'latex', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 22);
            
            % Panel E (Right - Tabs)
            rightPanel = uipanel(gridDisp, 'BorderType', 'none');
            tg = uitabgroup(rightPanel, 'Position', [0, 0, 300, 380]);
            
            % Tab 1
            tabStep = uitab(tg, 'Title', 'Last Step (E)');
            tabStep.BackgroundColor = [0 0 0]; 
            axE = uiaxes(tabStep, 'Position', [10, 10, 280, 320], 'Color', [0 0 0]);
            axis(axE, 'off');
            obj.MatrixDisplayE = text(axE, 0.5, 0.5, '', 'Interpreter', 'latex', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 16, 'Color', 'white');
            
            % Tab 2
            tabTotal = uitab(tg, 'Title', 'Total Product (P)');
            tabTotal.BackgroundColor = [0 0 0]; 
            axT = uiaxes(tabTotal, 'Position', [10, 10, 280, 320], 'Color', [0 0 0]);
            axis(axT, 'off');
            obj.MatrixDisplayTotal = text(axT, 0.5, 0.5, '', 'Interpreter', 'latex', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 14, 'Color', 'cyan'); 
            
            % 3. Command Center
            panelCmd = uipanel(obj.Fig, 'Title', 'Command Input', ...
                'Position', [20, 110, 850, 120]);
            
            gridCmd = uigridlayout(panelCmd, [2, 2]);
            gridCmd.RowHeight = {'1x', '1x'};
            gridCmd.ColumnWidth = {'4x', '1x'};
            
            obj.CmdInput = uieditfield(gridCmd, 'text', ...
                'Placeholder', 'Type command (e.g. R2 = R2 - R1) ... Press UP for History', ...
                'FontSize', 16, ...
                'ValueChangedFcn', @(txt,event) obj.processCommand()); 
            
            uibutton(gridCmd, 'Text', 'EXECUTE', ...
                'ButtonPushedFcn', @(btn,event) obj.processCommand(), ...
                'FontWeight', 'bold', 'BackgroundColor', [0.2 0.6 1], 'FontColor', 'white');
            
            lblHelp = uilabel(gridCmd, 'Text', 'Syntax: "R1 <-> R2", "R2 = R2 - 5*R1". (UP/DOWN Arrow to scroll history)');
            lblHelp.FontColor = [0.4 0.4 0.4];
            
            % 4. Substitution
            panelSub = uipanel(obj.Fig, 'Title', 'Substitute Variables', ...
                'Position', [20, 10, 850, 90]);
            gridSub = uigridlayout(panelSub, [1, 5]);
            
            uilabel(gridSub, 'Text', 'Replace:', 'HorizontalAlignment', 'right');
            obj.SubVarInput = uieditfield(gridSub, 'text', 'Placeholder', 'x');
            
            uilabel(gridSub, 'Text', 'With:', 'HorizontalAlignment', 'right');
            obj.SubValInput = uieditfield(gridSub, 'text', 'Placeholder', '5');
            
            uibutton(gridSub, 'Text', 'SUBSTITUTE', ...
                'ButtonPushedFcn', @(btn,event) obj.opSubstitute(), ...
                'BackgroundColor', [1 1 0.8]);

            % 5. History List
            panelLog = uipanel(obj.Fig, 'Title', 'History & Tools', ...
                'Position', [880, 10, 260, 620]);
            
            uibutton(panelLog, 'Text', 'Peek RREF', ...
                'Position', [10, 560, 240, 35], ...
                'BackgroundColor', [0.4 0.8 0.4], 'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(btn,event) obj.showRREF());
            
            obj.HistoryList = uilistbox(panelLog, 'Position', [10, 90, 240, 460]);
            
            uibutton(panelLog, 'Text', 'Undo', 'Position', [10, 50, 115, 30], ...
                'ButtonPushedFcn', @(btn,event) obj.opUndo());
            uibutton(panelLog, 'Text', 'Reset', 'Position', [135, 50, 115, 30], ...
                'ButtonPushedFcn', @(btn,event) obj.opReset());
        end
        
        % --- NEW: WINDOW KEY PRESS HANDLER ---
        function onWindowKeyPress(obj, event)
            % Check if the current focus is on our Command Input
            % (This prevents scrolling if you are focused on the listbox etc)
            
            if isempty(obj.CmdHistory), return; end
            
            % Only act if the focus is the EditField
            if obj.Fig.CurrentObject == obj.CmdInput
                
                switch event.Key
                    case 'uparrow'
                        % Move Index Backwards
                        obj.CmdHistoryIdx = max(1, obj.CmdHistoryIdx - 1);
                        if obj.CmdHistoryIdx <= length(obj.CmdHistory)
                             obj.CmdInput.Value = obj.CmdHistory(obj.CmdHistoryIdx);
                        end
                        
                    case 'downarrow'
                        % Move Index Forwards
                        obj.CmdHistoryIdx = min(length(obj.CmdHistory) + 1, obj.CmdHistoryIdx + 1);
                        
                        if obj.CmdHistoryIdx > length(obj.CmdHistory)
                            % If we go past the end, clear the box (like terminal)
                            obj.CmdInput.Value = "";
                        else
                            obj.CmdInput.Value = obj.CmdHistory(obj.CmdHistoryIdx);
                        end
                end
            end
        end

        % --- COMMAND PARSER ---
        function processCommand(obj)
            cmd = strtrim(string(obj.CmdInput.Value));
            if cmd == "", return; end
            
            % 1. Add to Command History (NEW)
            obj.CmdHistory(end+1) = cmd;
            obj.CmdHistoryIdx = length(obj.CmdHistory) + 1; % Reset pointer to end (blank)
            
            cmdUpper = upper(cmd);
            
            try
                if contains(cmdUpper, "<->") || contains(cmdUpper, "<=>")
                    parts = split(cmdUpper, ["<->", "<=>"]);
                    r1_str = regexp(parts(1), '\d+', 'match', 'once');
                    r2_str = regexp(parts(2), '\d+', 'match', 'once');
                    
                    if isempty(r1_str) || isempty(r2_str), error('Invalid Swap Syntax.'); end
                    
                    obj.performSwap(str2double(r1_str), str2double(r2_str));
                    
                elseif contains(cmdUpper, "=")
                    parts = split(cmd, "="); 
                    lhs = strtrim(parts(1));
                    targetIdx_str = regexp(lhs, '\d+', 'match', 'once');
                    
                    if isempty(targetIdx_str), error('Invalid Target Row.'); end
                    
                    obj.performRowOp(str2double(targetIdx_str), strtrim(parts(2))); 
                else
                    error('Unknown command.');
                end
                obj.CmdInput.Value = ''; % Clear input after success
                
            catch ME
                uialert(obj.Fig, ME.message, 'Command Error');
                % Don't clear input on error so user can fix typo
            end
        end
        
        % --- EXECUTION: SWAP ---
        function performSwap(obj, r1, r2)
            n = size(obj.CurrentMatrix, 1);
            if r1 > n || r2 > n || r1 < 1 || r2 < 1, error('Index out of bounds.'); end
            
            obj.saveState(sprintf('R%d <-> R%d', r1, r2));
            
            E = eye(n, 'sym');
            E([r1, r2], :) = E([r2, r1], :);
            obj.LastE = E;
            obj.TotalE = simplify(E * obj.TotalE);
            obj.CurrentMatrix = E * obj.CurrentMatrix;
            obj.updateDisplay();
        end
        
        % --- EXECUTION: ROW OP ---
        function performRowOp(obj, targetIdx, exprStr)
            n = size(obj.CurrentMatrix, 1);
            if targetIdx > n || targetIdx < 1, error('Target out of bounds.'); end
            
            rowTokens = sym('SYS_ROW_', [1 n]); 
            cleanExpr = exprStr;
            for i = 1:n
                pat = sprintf('(?<=\\W|^)[rR]%d(?=\\W|$)', i);
                cleanExpr = regexprep(cleanExpr, pat, sprintf('SYS_ROW_%d', i));
            end
            
            try, mathExpr = str2sym(cleanExpr); catch, error('Syntax Error'); end
            
            newRowA = zeros(1, size(obj.CurrentMatrix, 2), 'sym');
            newRowE = zeros(1, n, 'sym');
            I = eye(n, 'sym');
            
            for i = 1:n
                c_i = diff(mathExpr, rowTokens(i)); 
                if any(has(c_i, rowTokens)), error('Non-linear op.'); end
                if c_i ~= 0
                    newRowA = newRowA + c_i * obj.CurrentMatrix(i, :);
                    newRowE = newRowE + c_i * I(i, :);
                end
            end
            
            if subs(mathExpr, rowTokens, zeros(1,n)) ~= 0, error('No constants.'); end
            
            obj.saveState(sprintf('R%d = %s', targetIdx, exprStr));
            
            E = eye(n, 'sym');
            E(targetIdx, :) = simplify(newRowE);
            obj.LastE = E;
            obj.TotalE = simplify(E * obj.TotalE);
            obj.CurrentMatrix(targetIdx, :) = simplify(newRowA);
            obj.updateDisplay();
        end
        
        % --- RREF PEEK ---
        function showRREF(obj)
            d = uifigure('Name', 'RREF View', 'Position', [obj.Fig.Position(1)+50, obj.Fig.Position(2)+50, 600, 400]);
            lbl = uilabel(d, 'Text', 'Computing...', 'Position', [20, 350, 400, 30], 'FontSize', 16, 'FontWeight', 'bold');
            drawnow; 
            try
                R = rref(obj.CurrentMatrix); 
                ax = uiaxes(d, 'Position', [20, 20, 560, 320]); axis(ax, 'off');
                ax.XLim=[0 1]; ax.YLim=[0 1];
                l_str = latex(R);
                if length(l_str)>10000
                    lbl.Text='Too large to render.'; uitextarea(d,'Position',[20,20,560,320],'Value',char(R)); delete(ax);
                else
                    text(ax,0.5,0.5,['$$ \left[ ' l_str ' \right] $$'], 'Interpreter','latex','HorizontalAlignment','center','FontSize',18);
                    lbl.Text='Reduced Row Echelon Form:';
                end
            catch ME, lbl.Text='Error'; uialert(d, ME.message, 'Error'); end
        end
        
        % --- SUBSTITUTION ---
        function opSubstitute(obj)
            try
                sVar = sym(obj.SubVarInput.Value);
                sVal = sym(obj.SubValInput.Value);
                if isempty(sVar), return; end
                obj.saveState(sprintf('Sub: %s -> %s', string(sVar), string(sVal)));
                obj.LastE = eye(size(obj.CurrentMatrix,1), 'sym'); 
                obj.CurrentMatrix = subs(obj.CurrentMatrix, sVar, sVal);
                obj.TotalE = subs(obj.TotalE, sVar, sVal);
                obj.updateDisplay();
            catch ME, uialert(obj.Fig, ME.message, 'Sub Error'); end
        end
        
        % --- UNDO / RESET / DISPLAY ---
        function saveState(obj, desc)
            st.M = obj.CurrentMatrix; st.E = obj.LastE; st.T = obj.TotalE;
            obj.UndoStack{end+1} = st; obj.History{end+1} = desc;
        end
        function opUndo(obj)
            if isempty(obj.UndoStack), return; end
            st = obj.UndoStack{end}; obj.CurrentMatrix = st.M; obj.LastE = st.E; obj.TotalE = st.T;
            obj.UndoStack(end) = []; obj.History(end) = []; obj.updateDisplay();
        end
        function opReset(obj)
            obj.CurrentMatrix = obj.OriginalMatrix; n=size(obj.CurrentMatrix,1);
            obj.LastE = eye(n,'sym'); obj.TotalE = eye(n,'sym');
            obj.UndoStack={}; obj.History={'--- Reset ---'}; obj.updateDisplay();
            obj.CmdHistory = string.empty; obj.CmdHistoryIdx = 1; % Reset command history too
        end
        function updateDisplay(obj)
            try, obj.MatrixDisplayA.String = ['$$ \left[ ' latex(obj.CurrentMatrix) ' \right] $$']; catch, obj.MatrixDisplayA.String='Err'; end
            try, if isequal(obj.LastE, eye(size(obj.CurrentMatrix,1),'sym')), s='$$ I $$'; else, s=['$$ \left[ ' latex(obj.LastE) ' \right] $$']; end
            obj.MatrixDisplayE.String = s; catch, obj.MatrixDisplayE.String='?'; end
            try, obj.MatrixDisplayTotal.String = ['$$ \left[ ' latex(obj.TotalE) ' \right] $$']; catch, obj.MatrixDisplayTotal.String='Too Large'; end
            obj.HistoryList.Items = obj.History; if ~isempty(obj.History), scroll(obj.HistoryList, 'bottom'); end
        end
    end
end