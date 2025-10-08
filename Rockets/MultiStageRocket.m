classdef MultiStageRocket
    %MULTISTAGEROCKET Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name string
        stages (:, 1) % Needs 2n-1 rockets (and aerotables) where n is number of stages
        lockouts (:, :) 
        staging (:, 1)
    end
    
    methods
        function obj = MultiStageRocket(Name, stages)
            %MULTISTAGEROCKET Construct an instance of this class
            %   Detailed explanation goes here
            obj.Name = Name;
            obj.stages = stages;
        end

        function save(obj)
            %SAVE Save rocket to 'Rockets' folder in a convenient way.
            
            multistagerocket = obj;
            
            path = string("Rockets\" + obj.Name + ".mat");

            save(path, "multistagerocket")
        end

        function [] = display_geo(obj)
            %DISP Display wireframe of rocket geometry
            %   Show a representation of the rocket to verify that it looks
            %   practical and valid.
            
            fin_set_number = numel(obj.stages(1).FN1_height);
            stage_number = (numel(obj.stages) + 1) / 2;

            figure
            NC = polyshape([0,     1,      1] *  obj.stages(1).NC_length, ...
                           [0, 1 / 2, -1 / 2] *  obj.stages(1).NC_diam); % Need to acutally use correct nosecone shape
            stage_indices = (1:stage_number) * 2 - [(1:(stage_number - 1)) * 0, 1];
            for sn = 1:stage_number
                s = stage_indices(sn);
                AF_offset_rockets = obj.stages(stage_indices(sn:stage_number));
                AF_offset = -obj.stages(s).AF1_length;
                for o = 1:numel(AF_offset_rockets)
                    AF_offset = AF_offset + AF_offset_rockets(o).AF1_length;
                end
                AF1(sn) = polyshape([0,     1,    1,    0] *  obj.stages(s).AF1_length +  obj.stages(end).NC_length + AF_offset, ...
                                [ obj.stages(s).NC_diam / 2,  obj.stages(s).AF1_diam / 2, - obj.stages(s).AF1_diam / 2, - obj.stages(s).NC_diam / 2]);
            end
            for f = 1:fin_set_number
                FN1_top(f) = polyshape([0,  obj.stages(1).FN1_height(f) * tand(obj.stages(1).FN1_sweep(f)),  obj.stages(1).FN1_height(f) * tand(obj.stages(1).FN1_sweep(f)) +  obj.stages(1).FN1_tip(f),  obj.stages(1).FN1_root(f)] +  obj.stages(1).FN1_XLE(f), ...
                                    [0, 1, 1, 0] *  obj.stages(1).FN1_height(f) +  obj.stages(1).AF1_diam / 2);
                FN1_bottom(f) = polyshape([0,  obj.stages(1).FN1_height(f) * tand(obj.stages(1).FN1_sweep(f)),  obj.stages(1).FN1_height(f) * tand(obj.stages(1).FN1_sweep(f)) +  obj.stages(1).FN1_tip(f),  obj.stages(1).FN1_root(f)] +  obj.stages(1).FN1_XLE(f), ...
                                       -([0, 1, 1, 0] *  obj.stages(1).FN1_height(f) +  obj.stages(1).AF1_diam / 2));
            end
            plot(NC); hold on
            for s = 1:stage_number
                plot(AF1(s)); hold on
            end
            for f = 1:fin_set_number
                plot(FN1_top(f)); hold on
                plot(FN1_bottom(f)); hold on
            end
            scatter(obj.stages(1).COM_table(:, 2), zeros([size(obj.stages(1).COM_table, 1), 1]), 15, "blue", "filled"); hold on
            scatter(obj.stages(1).XCG_dry, 0,30,"red", "filled"); hold off
            axis equal
            xlabel("Axial [m]")
            ylabel("Transverse [m]")
            title("Rocket Diagram for " + obj.Name)
            legend([strings(1, stage_number + fin_set_number * 2 + 1), "COM table", "XCG dry"])
        end

        function [aerotables] = create_aerotables(obj, thread_number, target_batch_size)
            for s = 1:numel(obj.stages)
                aerotables(s).Aerotable = obj.stages(s).create_aerotable(thread_number, target_batch_size, Name_tag = obj.Name, aerotable_parameter_dialog = true);
            end
        end
    end

    methods(Static)
        function [obj] = stack(Name, individual_stages)
            n = numel(individual_stages);
            for s = 1:(n - 1)
                stages(2 * s - 1) = individual_stages(s + 1) + individual_stages(s);
                stages(2 * s) = individual_stages(s);
                stages(2 * s) = stages(2 * s).modify(array2table([0], "VariableNames", "NC_length"));
            end
            stages(2 * n - 1) = individual_stages(end);

            obj = MultiStageRocket(Name, stages);
        end
    end
end

