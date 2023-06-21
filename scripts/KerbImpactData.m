classdef KerbImpactData
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties(Constant)
        DATA_START_IDX = 43;
        get_data_vector = @(raw_data, channel_idx) raw_data(channel_idx).data(KerbImpactData.DATA_START_IDX:end);
        get_data_point = @(raw_data , channel_idx , point_idx) raw_data(channel_idx).data(point_idx);
    end
    properties
        run_name
        run_speed_kph
        run_tyre_pressure_bar
        plot_number
        asc_file_name
        raw_time
        shifted_time
        drawstring_position
        drawstring_speed
        drawstring_speed_kph
        run_date_time
        trigger
        trigger_idx
        brake_pressure
        fl
        fr
        rr
        rl
        end_idx
    end

    methods
        function obj = KerbImpactData(data_cell,row_num, channel_idx)
            % makes a data object based on the carData cell array 
            % indices for channels are hard coded,
            % shifted time is zero at trigger time
            raw_data = {data_cell{row_num , :}};
            obj.run_name = raw_data{2};
            obj.asc_file_name = raw_data{4};
            obj = obj.parse_data(raw_data , channel_idx);
        end
        function obj = parse_data(obj, raw_data , channel_idx)
            obj.raw_time = obj.get_data_vector(raw_data{1} , channel_idx.time);
            obj.trigger = obj.get_data_vector(raw_data{1}, channel_idx.camTrigger);
            obj.trigger_idx = find(obj.trigger>0.5, 1 ,'first');
            obj = obj.parse_config();
            obj = obj.parse_logger_data(raw_data{1}, channel_idx);
            obj = obj.parse_drawstring_data(raw_data );
            obj = obj.set_end_idx();
            %%% last step to shift position to be measured from the step
            %%% point in laser
            [~ , step_idx] = max(obj.fl.laser_height); 
            obj.drawstring_position = obj.drawstring_position - obj.drawstring_position(step_idx);
        end
        function obj = parse_logger_data(obj , raw_data , channel_idx)
            obj.shifted_time = obj.raw_time - obj.raw_time(obj.trigger_idx);
            obj.run_date_time = obj.get_time_date(raw_data, channel_idx, obj.trigger_idx);
            obj.brake_pressure = obj.get_data_vector(raw_data , channel_idx.brake_pressure);
            obj.fl = obj.get_corner_data(raw_data, channel_idx, 'fl');
            obj.fr = obj.get_corner_data(raw_data, channel_idx, 'fr');
            obj.rl = obj.get_corner_data(raw_data, channel_idx, 'rl');
            obj.rr = obj.get_corner_data(raw_data, channel_idx, 'rr');
        end
        function obj = parse_drawstring_data(obj, raw_data)
            draw_string_time = raw_data{3}(: , 1);
            draw_string_position = raw_data{3}(: , 2);
            draw_string_trigger = raw_data{3}(: , 3);
            draw_string_trigger_idx = find(draw_string_trigger>0.5 , 1 , 'first');
            draw_string_time = draw_string_time - draw_string_time(draw_string_trigger_idx);
            obj.drawstring_position = interp1(draw_string_time, draw_string_position, obj.shifted_time);
            obj.drawstring_speed = diff(obj.drawstring_position)./diff(obj.shifted_time);
            % keeping the length the same
            obj.drawstring_speed(end+1) = obj.drawstring_speed(end);
            obj.drawstring_speed_kph = obj.drawstring_speed * 3.6;
        end
        function obj = parse_config(obj)
            pressure_string = regexp(obj.run_name, "\dp\d", 'match');
            obj.run_tyre_pressure_bar = str2double(pressure_string{1}(1)) + str2double(pressure_string{1}(3))*0.1;
            speed_string = regexp(obj.run_name,"\d*kph", 'match');
            obj.run_speed_kph = str2double(regexp(speed_string{1} , "\d*", 'match'));
            obj.plot_number = 1 + (obj.run_speed_kph == 7)*2 + (obj.run_tyre_pressure_bar == 2.6);
        end
        function obj = set_end_idx(obj)
            start_idx = find(obj.drawstring_speed_kph > (obj.run_speed_kph-2) , 1 , "first");
            obj.end_idx = find(obj.brake_pressure(start_idx:end)> 10 , 1 , 'first') + start_idx;
        end
        function line_obj = quad_plot(obj , corner_name, value_name, line_spec)
            if nargin <4
                line_spec = "-b";
            end
            subplot(4 , 1 , obj.plot_number);
            hold on 
            line_obj = plot(obj.drawstring_position(1:obj.end_idx), obj.(corner_name).(value_name)(1:obj.end_idx), line_spec);
            title (sprintf("pressure: %.1f bar, speed: %.0f kph", obj.run_tyre_pressure_bar, obj.run_speed_kph));
        end
    end
        methods(Static)
            function date_time = get_time_date(raw_data, channel_idx, trigger_idx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            date_time = datetime(...
                KerbImpactData.get_data_point(raw_data, channel_idx.year , trigger_idx),...
                KerbImpactData.get_data_point(raw_data, channel_idx.month , trigger_idx),...
                KerbImpactData.get_data_point(raw_data, channel_idx.day , trigger_idx),...
                KerbImpactData.get_data_point(raw_data, channel_idx.hour , trigger_idx),...
                KerbImpactData.get_data_point(raw_data, channel_idx.minute , trigger_idx),...
                KerbImpactData.get_data_point(raw_data, channel_idx.second, trigger_idx));
            end
            function corner_data = get_corner_data(raw_data, channel_idx, corner_name)
                signal_names = fieldnames(channel_idx.(corner_name));
                for i = 1:length(signal_names)
                    name = signal_names{i};
                    corner_data.(name) = KerbImpactData.get_data_vector(raw_data , channel_idx.(corner_name).(name));
                end
            end
end
end