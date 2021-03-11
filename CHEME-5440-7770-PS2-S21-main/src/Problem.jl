function generate_problem_dictionary(path_to_parameters_file::String)::Dict{String,Any}

    # initialize -
    problem_dictionary = Dict{String,Any}()

    try

        # load the TOML parameters file -
        toml_dictionary = TOML.parsefile(path_to_parameters_file)["biophysical_constants"]

        # setup the initial condition array -
        initital_condition_array = [
            0.0 ;   # 1 mRNA
            .005 ;   # TODO: gene concentration goes here -
            0.0 ;   # 3 I = we'll fill this in the execute script 
        ]


        # TODO: calculate the mRNA_degradation_constant 
        # mRNA_degradation_constant = ....

        mRNA_degradation_constant = log(2)/toml_dictionary["mRNA_half_life_in_min"]

        # TODO: VMAX for transcription -
        # VMAX = ...

        VMAX = (toml_dictionary["transcription_elongation_rate"]/toml_dictionary["gene_length_in_nt"])*toml_dictionary["RNAPII_concentration"]

        # TODO: Stuff that I'm forgetting?

        # --- PUT STUFF INTO problem_dictionary ---- 
        problem_dictionary["transcription_time_constant"] = toml_dictionary["transcription_time_constant"]
        problem_dictionary["transcription_saturation_constant"] = toml_dictionary["transcription_saturation_constant"]
        problem_dictionary["E1"] = toml_dictionary["energy_promoter_state_1"]
        problem_dictionary["E2"] = toml_dictionary["energy_promoter_state_2"]
        problem_dictionary["inducer_dissociation_constant"] = toml_dictionary["inducer_dissociation_constant"]
        problem_dictionary["inducer_cooperativity_parameter"] = toml_dictionary["inducer_cooperativity_parameter"]
        problem_dictionary["ideal_gas_constant_R"] = 0.008314 # kJ/mol-K
        problem_dictionary["temperature_K"] = (273.15+37)
        problem_dictionary["initial_condition_array"] = initital_condition_array
        problem_dictionary["maximum_transcription_velocity"] = VMAX
        problem_dictionary["mRNA_degradation_constant"] = mRNA_degradation_constant
        
        # return -
        return problem_dictionary
    catch error
        throw(error)
    end
end