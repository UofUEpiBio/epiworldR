#ifndef EPIWORLD_MODEL_MEAT_PRINT_HPP
#define EPIWORLD_MODEL_MEAT_PRINT_HPP

template<typename TSeq>
inline const Model<TSeq> & Model<TSeq>::print(bool lite) const
{

    // Horizontal line
    std::string line = "";
    for (epiworld_fast_uint i = 0u; i < 80u; ++i)
        line += "_";

    // Prints a message if debugging is on
    EPI_DEBUG_NOTIFY_ACTIVE()

    printf_epiworld("%s\n",line.c_str());

    if (lite)
    {
        // Printing the name of the model
        printf_epiworld("%s", name.c_str());

        // Printing the number of agents, viruses, and tools
        printf_epiworld(
            "\nIt features %i agents, %i virus(es), and %i tool(s).\n",
            static_cast<int>(size()),
            static_cast<int>(get_n_viruses()),
            static_cast<int>(get_n_tools())
            );

        printf_epiworld(
            "The model has %i states.",
            static_cast<int>(nstates)
            );

        if (today() != 0)
        {
            printf_epiworld(
                "\nThe final distribution is: "
            );

            int nstate_int = static_cast<int>(nstates);

            for (int i = 0u; i < nstate_int; ++i)
            {
                printf_epiworld(
                    "%i %s%s",
                    static_cast<int>(db.today_total[ i ]),
                    states_labels[i].c_str(),
                    (
                        i == (nstate_int - 2)
                        ) ? ", and " : (
                            (i == (nstate_int - 1)) ? ".\n" : ", "
                            )
                );
            }
        } else {
            printf_epiworld(" The model hasn't been run yet.\n");
        }

        return *this;
    }

    printf_epiworld("%s\n%s\n\n",line.c_str(), "SIMULATION STUDY");

    printf_epiworld("Name of the model   : %s\n", (this->name == "") ? std::string("(none)").c_str() : name.c_str());
    printf_epiworld("Population size     : %i\n", static_cast<int>(size()));

    auto ncols = get_agents_data_ncols();

    if (ncols > 0)
    {
        printf_epiworld("Agents' data loaded : yes (%i columns/features)\n", static_cast<int>(ncols));
    }
    else
    {
        printf_epiworld("Agents' data        : (none)\n");
    }

    printf_epiworld("Number of entities  : %i\n", static_cast<int>(entities.size()));
    printf_epiworld("Days (duration)     : %i (of %i)\n", today(), static_cast<int>(ndays));
    printf_epiworld("Number of viruses   : %i\n", static_cast<int>(db.get_n_viruses()));
    if (n_replicates > 0u)
    {
        std::string abbr;
        epiworld_double elapsed;
        epiworld_double total;
        get_elapsed("auto", &elapsed, &total, &abbr, false);
        printf_epiworld("Last run elapsed t  : %.2f%s\n", elapsed, abbr.c_str());
        if (n_replicates > 1u)
        {
            printf_epiworld("Total elapsed t     : %.2f%s (%i runs)\n", total, abbr.c_str(), static_cast<int>(n_replicates));
        }

        // Elapsed time in speed
        get_elapsed("microseconds", &elapsed, &total, &abbr, false);
        printf_epiworld("Last run speed      : %.2f million agents x day / second\n",
            static_cast<double>(this->size()) *
            static_cast<double>(this->get_ndays()) /
            static_cast<double>(elapsed)
            );
        if (n_replicates > 1u)
        {
            printf_epiworld("Average run speed   : %.2f million agents x day / second\n",
                static_cast<double>(this->size()) *
                static_cast<double>(this->get_ndays()) *
                static_cast<double>(n_replicates) /
                static_cast<double>(total)
            );
        }

    } else {
        printf_epiworld("Last run elapsed t  : -\n");
    }
    
    
    if (rewire_fun)
    {
        printf_epiworld("Rewiring            : on (%.2f)\n\n", rewire_prop);
    } else {
        printf_epiworld("Rewiring            : off\n\n");
    }
    
    // Printing Global events
    printf_epiworld("Global events:\n");
    for (auto & a : globalevents)
    {
        if (a.get_day() < 0)
        {
            printf_epiworld(" - %s (runs daily)\n", a.get_name().c_str());
        } else {
            printf_epiworld(" - %s (day %i)\n", a.get_name().c_str(), a.get_day());
        }
    }

    if (globalevents.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    printf_epiworld("\nVirus(es):\n");
    size_t n_viruses_model = viruses.size();
    for (size_t i = 0u; i < n_viruses_model; ++i)
    {    

        
        const auto & virus = viruses[i];
        if ((n_viruses_model > 10) && (i >= 10))
        {
            printf_epiworld(" ...and %i more viruses...\n",
                static_cast<int>(n_viruses_model) - 
                static_cast<int>(i)
                );
            break;
        }

        if (i < n_viruses_model)
        {

            printf_epiworld(
                " - %s\n",
                virus->get_name().c_str()
            );

        } else {

            printf_epiworld(
                " - %s (originated in the model...)\n",
                virus->get_name().c_str()
            );

        }

    }

    auto nvariants = db.get_n_viruses() - n_viruses_model;
    if (nvariants > 0)
    {

        printf_epiworld(" ...and %i more variants...\n", static_cast<int>(nvariants));

    }

    if (viruses.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    printf_epiworld("\nTool(s):\n");
    size_t n_tools_model = tools.size();
    for (size_t i = 0u; i < tools.size(); ++i)
    {   
        const auto & tool = tools[i];

        if ((n_tools_model > 10) && (i >= 10))
        {
            printf_epiworld(
                " ...and %i more tools...\n",
                static_cast<int>(n_tools_model) - static_cast<int>(i)
                );
            break;
        }

        if (i < n_tools_model)
        {
            printf_epiworld(
                " - %s\n",
                tool->get_name().c_str()
                );


        } else {

            printf_epiworld(
                " - %s (originated in the model...)\n",
                tool->get_name().c_str()
            );

        }
        

    }

    if (tools.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    // Information about the parameters included
    printf_epiworld("\nModel parameters:\n");
    epiworld_fast_uint nchar = 0u;
    for (auto & p : parameters)
        if (p.first.length() > nchar)
            nchar = p.first.length();

    std::string fmt = " - %-" + std::to_string(nchar + 1) + "s: ";
    for (auto & p : parameters)
    {
        std::string fmt_tmp = fmt;
        if (std::fabs(p.second) < 0.0001)
            fmt_tmp += "%.1e\n";
        else
            fmt_tmp += "%.4f\n";

        printf_epiworld(
            fmt_tmp.c_str(),
            p.first.c_str(),
            p.second
        );
        
    }

    if (parameters.size() == 0u)
    {
        printf_epiworld(" (none)\n");
    }

    nchar = 0u;
    for (auto & p : states_labels)
        if (p.length() > nchar)
            nchar = p.length();

    

    if (today() != 0)
    {
        fmt =
            std::string("  - (%") +
                std::to_string(std::to_string(nstates).length()) +
            std::string("d) %-") + std::to_string(nchar) +
            std::string("s : %") +
            std::to_string(std::to_string(size()).length()) +
            std::string("i -> %i\n");
    } else {
        fmt =
            std::string("  - (%") +
                std::to_string(std::to_string(nstates).length()) +
            std::string("d) %-") + std::to_string(nchar) +
            std::string("s : %i\n");
    }
        
    if (today() != 0)
    {
        printf_epiworld("\nDistribution of the population at time %i:\n", today());
        for (size_t s = 0u; s < nstates; ++s)
        {

                printf_epiworld(
                    fmt.c_str(),
                    s,
                    states_labels[s].c_str(),
                    db.hist_total_counts[s],
                    db.today_total[ s ]
                    );

        }
    }

    if (today() != 0)
        (void) db.get_transition_probability(true);

    return *this;

}

#endif