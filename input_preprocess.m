function [input, output, mean, std, target, test_stations] = input_preprocess(Data, m, s, PM10,...
    NumOutputs, CoordFlag, Tval, Nstaz_val, ValSelection)

switch ValSelection
    case 0
        % N_input = 3;
        
        % Prepare input
        % for i = 1:N_input-2
        %     Data_id{i} = Data{i}(1:end-NumOutputs,:);
        % end
        
        Data_id{1} = Data{1}(1:end-NumOutputs,:);
        
        % Prepare target
        for i = 1:NumOutputs
            Data_tar{i} = Data{1}(1+i:end-NumOutputs+i, :);
            m_tar{i} = m{1}(1+i:end-NumOutputs+i, :);
            s_tar{i} = s{1}(1+i:end-NumOutputs+i, :);
            PM10_target{i} = PM10(1+i:end-NumOutputs+i, :);
        end
        
        %% Identification dataset
        Nstaz = width(Data{1,1})-Nstaz_val;
        N = (length(Data_id{1})-Tval)*Nstaz;
        
        idx_row_id = false(length(Data_id{1}), 1);
        idx_row_id(1:end-Tval) = 1;
        idx_col_id = true(width(Data_id{1}), 1);
        idx = randperm(width(Data_id{1}), width(Data_id{1})-Nstaz);
        test_stations = idx;
        idx_col_id(idx) = 0;
        
        % for i = 1:N_input-2
        %     input_id{i} = reshape(Data_id{i}(idx_row_id,idx_col_id), N, 1);
        % end
        input_id{1} = reshape(Data_id{1}(idx_row_id,idx_col_id), N, 1);
        
        
        if CoordFlag == true
            input_id{2} = repmat(Data{end-1}(idx_col_id), (length(Data_id{1})-Tval), 1);
            input_id{3} = repmat(Data{end}(idx_col_id), (length(Data_id{1})-Tval), 1);
            input_id{2} = reshape(input_id{2}, N, 1);
            input_id{3} = reshape(input_id{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_id{i} = reshape(Data_tar{i}(idx_row_id,idx_col_id), N, 1);
            output_id{i} = output_id{i}';
            m_tar_id{i} = reshape(m_tar{i}(idx_row_id,idx_col_id), N, 1);
            m_tar_id{i} = m_tar_id{i}';
            s_tar_id{i} = reshape(s_tar{i}(idx_row_id,idx_col_id), N, 1);
            s_tar_id{i} = s_tar_id{i}';
            PM10_target_id{i} = reshape(PM10_target{i}(idx_row_id,idx_col_id), N, 1);
            PM10_target_id{i} = PM10_target_id{i}';
        end
        
        %% Validation dataset
        N = (Tval)*Nstaz;
        
        idx_row_val = true(length(Data_id{1}), 1);
        idx_row_val(1:end-Tval) = 0;
        
        % for i = 1:N_input-2
        %     input_Tval{i} = reshape(Data_id{i}(idx_row_Tval,idx_col_id), N, 1);
        % end
        input_val{1} = reshape(Data_id{1}(idx_row_val,idx_col_id), N, 1);
        
        
        if CoordFlag == true
            input_val{2} = repmat(Data{end-1}(idx_col_id), (Tval), 1);
            input_val{3} = repmat(Data{end}(idx_col_id), (Tval), 1);
            input_val{2} = reshape(input_val{2}, N, 1);
            input_val{3} = reshape(input_val{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_val{i} = reshape(Data_tar{i}(idx_row_val,idx_col_id), N, 1);
            output_val{i} = output_val{i}';
            m_tar_val{i} = reshape(m_tar{i}(idx_row_val,idx_col_id), N, 1);
            m_tar_val{i} = m_tar_val{i}';
            s_tar_val{i} = reshape(s_tar{i}(idx_row_val,idx_col_id), N, 1);
            s_tar_val{i} = s_tar_val{i}';
            PM10_target_val{i} = reshape(PM10_target{i}(idx_row_val,idx_col_id), N, 1);
            PM10_target_val{i} = PM10_target_val{i}';
        end
        
        %% Test dataset
        N = length(Data_id{1})*length(idx);
        idx_col_test = false(width(Data_id{1}), 1);
        idx_col_test(idx) = 1;
        
        % for i = 1:N_input-2
        %     input_Sval{i} = reshape(Data_id{i}(:,idx_col_val), N, 1);
        % end
        input_test{1} = reshape(Data_id{1}(:,idx_col_test), N, 1);
        
        if CoordFlag == true
            input_test{2} = repmat(Data{end-1}(idx_col_test), (length(Data_id{1})), 1);
            input_test{3} = repmat(Data{end}(idx_col_test), (length(Data_id{1})), 1);
            input_test{2} = reshape(input_test{2}, N, 1);
            input_test{3} = reshape(input_test{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_test{i} = reshape(Data_tar{i}(:,idx_col_test), N, 1);
            output_test{i} = output_test{i}';
            m_tar_test{i} = reshape(m_tar{i}(:,idx_col_test), N, 1);
            m_tar_test{i} = m_tar_test{i}';
            s_tar_test{i} = reshape(s_tar{i}(:,idx_col_test), N, 1);
            s_tar_test{i} = s_tar_test{i}';
            PM10_target_test{i} = reshape(PM10_target{i}(:,idx_col_test), N, 1);
            PM10_target_test{i} = PM10_target_test{i}';
        end
        
        %% Organize output
        input{1} = input_id;
        input{2} = input_val;
        input{3} = input_test;
        
        output{1} = output_id;
        output{2} = output_val;
        output{3} = output_test;
        
        mean{1} = m_tar_id;
        mean{2} = m_tar_val;
        mean{3} = m_tar_test;
        
        std{1} = s_tar_id;
        std{2} = s_tar_val;
        std{3} = s_tar_test;
        
        target{1} = PM10_target_id;
        target{2} = PM10_target_val;
        target{3} = PM10_target_test;

    case 1

        Tval = Tval*2;
        Data_id{1} = Data{1}(1:end-NumOutputs,:);
        
        % Prepare target
        for i = 1:NumOutputs
            Data_tar{i} = Data{1}(1+i:end-NumOutputs+i, :);
            m_tar{i} = m{1}(1+i:end-NumOutputs+i, :);
            s_tar{i} = s{1}(1+i:end-NumOutputs+i, :);
            PM10_target{i} = PM10(1+i:end-NumOutputs+i, :);
        end
        
        %% Input data and target
        Nstaz = width(Data{1,1});
        N = (length(Data_id{1})-Tval)*Nstaz;
        
        idx_row_id = false(length(Data_id{1}), 1);
        idx_row_id(1:end-Tval) = 1;
        test_stations = [];
        
        input_id{1} = reshape(Data_id{1}(idx_row_id,:), N, 1);
        
        
        if CoordFlag == true
            % input_id{2} = repmat(Data{end-1}(:), (length(Data_id{1})-Tval), 1);
            % input_id{3} = repmat(Data{end}(:), (length(Data_id{1})-Tval), 1);
            % input_id{2} = reshape(input_id{2}, N, 1);
            % input_id{3} = reshape(input_id{3}, N, 1);

            input_id{2} = repmat(Data{end-1}(:), (length(Data_id{1})-Tval), 1);
            input_id{3} = repmat(Data{end}(:), (length(Data_id{1})-Tval), 1);
            input_id{2} = reshape(input_id{2}, Nstaz, length(Data_id{1})-Tval);
            input_id{2} = input_id{2}';
            input_id{2} = reshape(input_id{2}, N, 1);
            input_id{3} = reshape(input_id{3}, Nstaz, length(Data_id{1})-Tval);
            input_id{3} = input_id{3}';
            input_id{3} = reshape(input_id{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_id{i} = reshape(Data_tar{i}(idx_row_id,:), N, 1);
            output_id{i} = output_id{i}';
            m_tar_id{i} = reshape(m_tar{i}(idx_row_id,:), N, 1);
            m_tar_id{i} = m_tar_id{i}';
            s_tar_id{i} = reshape(s_tar{i}(idx_row_id,:), N, 1);
            s_tar_id{i} = s_tar_id{i}';
            PM10_target_id{i} = reshape(PM10_target{i}(idx_row_id,:), N, 1);
            PM10_target_id{i} = PM10_target_id{i}';
        end

        %% Validation
        Tval = Tval/2;

        N = (Tval)*Nstaz;
        
        idx_row_val = false(length(Data_id{1}), 1);
        idx_row_val(end-(Tval*2)+1:end-Tval) = 1;
        
        % for i = 1:N_input-2
        %     input_Tval{i} = reshape(Data_id{i}(idx_row_Tval,idx_col_id), N, 1);
        % end
        input_val{1} = reshape(Data_id{1}(idx_row_val,:), N, 1);
        
        
        if CoordFlag == true
            % input_val{2} = repmat(Data{end-1}(:), (Tval), 1);
            % input_val{3} = repmat(Data{end}(:), (Tval), 1);
            % input_val{2} = reshape(input_val{2}, N, 1);
            % input_val{3} = reshape(input_val{3}, N, 1);

            input_val{2} = repmat(Data{end-1}(:), (Tval), 1);
            input_val{3} = repmat(Data{end}(:), (Tval), 1);
            input_val{2} = reshape(input_val{2}, Nstaz, Tval);
            input_val{2} = input_val{2}';
            input_val{2} = reshape(input_val{2}, N, 1);
            input_val{3} = reshape(input_val{3}, Nstaz, Tval);
            input_val{3} = input_val{3}';
            input_val{3} = reshape(input_val{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_val{i} = reshape(Data_tar{i}(idx_row_val,:), N, 1);
            output_val{i} = output_val{i}';
            m_tar_val{i} = reshape(m_tar{i}(idx_row_val,:), N, 1);
            m_tar_val{i} = m_tar_val{i}';
            s_tar_val{i} = reshape(s_tar{i}(idx_row_val,:), N, 1);
            s_tar_val{i} = s_tar_val{i}';
            PM10_target_val{i} = reshape(PM10_target{i}(idx_row_val,:), N, 1);
            PM10_target_val{i} = PM10_target_val{i}';
        end

        %% Test dataset
        N = (Tval)*Nstaz;
        idx_row_test = true(length(Data_id{1}), 1);
        idx_row_test(1:end-Tval) = 0;
        
        input_test{1} = reshape(Data_id{1}(idx_row_test,:), N, 1);
        
        
        if CoordFlag == true
            % input_test{2} = repmat(Data{end-1}(:), (Tval), 1);
            % input_test{3} = repmat(Data{end}(:), (Tval), 1);
            % input_test{2} = reshape(input_test{2}, N, 1);
            % input_test{3} = reshape(input_test{3}, N, 1);

            input_test{2} = repmat(Data{end-1}(:), (Tval), 1);
            input_test{3} = repmat(Data{end}(:), (Tval), 1);
            input_test{2} = reshape(input_test{2}, Nstaz, Tval);
            input_test{2} = input_test{2}';
            input_test{2} = reshape(input_test{2}, N, 1);
            input_test{3} = reshape(input_test{3}, Nstaz, Tval);
            input_test{3} = input_test{3}';
            input_test{3} = reshape(input_test{3}, N, 1);
        end
        
        for i = 1:NumOutputs
            output_test{i} = reshape(Data_tar{i}(idx_row_test,:), N, 1);
            output_test{i} = output_test{i}';
            m_tar_test{i} = reshape(m_tar{i}(idx_row_test,:), N, 1);
            m_tar_test{i} = m_tar_test{i}';
            s_tar_test{i} = reshape(s_tar{i}(idx_row_test,:), N, 1);
            s_tar_test{i} = s_tar_test{i}';
            PM10_target_test{i} = reshape(PM10_target{i}(idx_row_test,:), N, 1);
            PM10_target_test{i} = PM10_target_test{i}';
        end

        %% Organize output
        input{1} = input_id;
        input{2} = input_val;
        input{3} = input_test;
        
        output{1} = output_id;
        output{2} = output_val;
        output{3} = output_test;
        
        mean{1} = m_tar_id;
        mean{2} = m_tar_val;
        mean{3} = m_tar_test;
        
        std{1} = s_tar_id;
        std{2} = s_tar_val;
        std{3} = s_tar_test;
        
        target{1} = PM10_target_id;
        target{2} = PM10_target_val;
        target{3} = PM10_target_test;

end
