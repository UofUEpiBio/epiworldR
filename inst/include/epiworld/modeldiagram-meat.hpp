#ifndef EPIWORLD_MODELDIAGRAM_MEAT_HPP
#define EPIWORLD_MODELDIAGRAM_MEAT_HPP

inline void ModelDiagram::read_transitions(
    const std::string & fn_transition
) {

    // Checking if the file exists
    std::ifstream file(fn_transition);

    if (!file.is_open())
        throw std::runtime_error(
            "Could not open the file '" +
            fn_transition + "' for reading."
        );

    // Reading the data
    std::string line;
    int i = 0;
    while (std::getline(file, line))
    {

        // Skipping the first line
        if (i++ == 0)
            continue;

        std::istringstream iss(line);
        #ifdef EPI_DEBUG
        int t;
        iss >> t;
        #endif

        int i_;
        std::string from_, to_;
        int counts_;

        iss >> i_;

        // Read the quoted strings
        iss >> std::quoted(from_) >> std::quoted(to_);

        // Read the integer
        iss >> counts_;

        if (counts_ > 0)
        {
            auto idx = std::make_pair(from_, to_);
            if (data.find(idx) == data.end())
                data[idx] = counts_;
            else
                data[idx] += counts_;
        }

    }

    // Incrementing the number of runs
    this->n_runs++;

}

inline void ModelDiagram::read_transitions(
    const std::vector< std::string > & fns_transition
)
{

    for (const auto & fn: fns_transition)
        this->read_transitions(fn);

    return;

}

inline void ModelDiagram::transition_probability(
    bool normalize
)
{

    // Generating the map of states
    std::set< std::string > states_set;

    for (const auto & d: data)
    {
        states_set.insert(d.first.first);
        states_set.insert(d.first.second);
    }

    // Generating the transition matrix
    states = std::vector< std::string >(states_set.begin(), states_set.end());
    size_t n_states = states.size();
    tprob.resize(n_states * n_states);
    std::fill(tprob.begin(), tprob.end(), 0.0);

    std::vector< epiworld_double > rowsum(n_states, 0.0);
    for (size_t i = 0; i < n_states; ++i)
    {

        for (size_t j = 0; j < n_states; ++j)
        {

            auto key = std::make_pair(states[i], states[j]);
            if (data.find(key) != data.end())
                tprob[i + j * n_states] = static_cast<epiworld_double>(
                    data[key]
                );

            if (normalize)
                rowsum[i] += tprob[i + j * n_states];

        }

        if (normalize)
        {
            for (size_t j = 0; j < n_states; ++j)
                tprob[i + j * n_states] /= rowsum[i];
        }

    }

    return;


}

inline void ModelDiagram::clear()
{
    data.clear();
    states.clear();
    tprob.clear();
    n_runs = 0;
    return;
}

inline void ModelDiagram::draw_mermaid(
    std::string fn_output,
    bool self
) {
     // Getting a sorting vector of indices from the states
    // string vector
    std::vector< size_t > idx(states.size());
    std::iota(idx.begin(), idx.end(), 0u);

    std::sort(
        idx.begin(),
        idx.end(),
        [&states = this->states](size_t i, size_t j) {
            return states[i] < states[j];
        }
    );

    std::vector< std::string > states_ids;
    for (size_t i = 0u; i < states.size(); ++i)
        states_ids.push_back("s" + std::to_string(i));

    std::string graph = "flowchart LR\n";

    // Declaring the states
    for (size_t i = 0u; i < states.size(); ++i)
    {
        graph += "\t" + states_ids[i] + "[" + states[idx[i]] + "]\n";
    }

    // Adding the transitions
    size_t n_states = states.size();
    for (size_t i = 0u; i < states.size(); ++i)
    {
        for (size_t j = 0u; j < states.size(); ++j)
        {
            if (!self && i == j)
                continue;

            if (tprob[idx[i] + idx[j] * n_states] > 0.0)
            {
                graph += "\t" + states_ids[i] + " -->|" +
                    std::to_string(tprob[idx[i] + idx[j] * n_states]) +
                    "| " + states_ids[j] + "\n";
            }
        }
    }

    if (fn_output != "")
    {
        std::ofstream file(fn_output);

        if (!file.is_open())
            throw std::runtime_error(
                "Could not open the file " +
                fn_output +
                " for writing."
            );

        file << graph;
        file.close();

    } else {
        printf_epiworld("%s\n", graph.c_str());
    }

    return;
}

inline void ModelDiagram::draw_dot(std::string fn_output, bool self) {
    std::vector<size_t> idx(states.size());
    std::iota(idx.begin(), idx.end(), 0u);

    std::sort(
        idx.begin(),
        idx.end(),
        [&states = this->states](size_t i, size_t j) {
            return states[i] < states[j];
        }
    );

    std::vector<std::string> states_ids;
    for (size_t i = 0u; i < states.size(); ++i)
        states_ids.push_back("s" + std::to_string(i));

    std::string graph = "digraph G {\n\trankdir=LR;\n";

    // Declaring the states
    for (size_t i = 0u; i < states.size(); ++i)
    {
        graph += "\t" + states_ids[i] + " [label=\"" + states[idx[i]] + "\"];\n";
    }

    // Adding the transitions
    size_t n_states = states.size();
    for (size_t i = 0u; i < states.size(); ++i)
    {
        for (size_t j = 0u; j < states.size(); ++j)
        {
            if (!self && i == j)
                continue;

            if (tprob[idx[i] + idx[j] * n_states] > 0.0)
            {
                graph += "\t" + states_ids[i] + " -> " + states_ids[j] +
                         " [label=\"" +
                         std::to_string(tprob[idx[i] + idx[j] * n_states]) +
                         "\"];\n";
            }
        }
    }

    graph += "}\n";

    if (fn_output != "")
    {
        std::ofstream file(fn_output);

        if (!file.is_open())
            throw std::runtime_error(
                "Could not open the file " +
                fn_output +
                " for writing."
            );

        file << graph;
        file.close();

    } else {
        printf_epiworld("%s\n", graph.c_str());
    }

    return;
}

inline void ModelDiagram::draw(
    DiagramType diagram_type,
    std::string fn_output,
    bool self
)
{
    switch (diagram_type) {
    case DiagramType::Mermaid:
        return draw_mermaid(fn_output, self);
    case DiagramType::DOT:
        return draw_dot(fn_output, self);
    default:
  		throw std::logic_error("Unknown how to dispatch DiagramType.");
    }
}

inline void ModelDiagram::draw_from_file(
    DiagramType diagram_type,
    const std::string & fn_transition,
    const std::string & fn_output,
    bool self
) {

    this->clear();

    // Loading the transition file
    this->read_transitions(fn_transition);

    // Computing the transition probability
    this->transition_probability();

    // Actually drawing the diagram
    this->draw(diagram_type, fn_output, self);

    return;

}

inline void ModelDiagram::draw_from_files(
    DiagramType diagram_type,
    const std::vector< std::string > & fns_transition,
    const std::string & fn_output,
    bool self
) {

    this->clear();

    // Loading the transition files
    this->read_transitions(fns_transition);

    // Computing the transition probability
    this->transition_probability();

    // Actually drawing the diagram
    this->draw(diagram_type, fn_output, self);

    return;
}

inline void ModelDiagram::draw_from_data(
    DiagramType diagram_type,
    const std::vector< std::string > & states,
    const std::vector< epiworld_double > & tprob,
    const std::string & fn_output,
    bool self
) {

    this->clear();

    this->states = states;
    this->tprob = tprob;

    this->draw(diagram_type, fn_output, self);

    return;

}

#endif
