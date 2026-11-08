#ifndef EPIWORLD_REGISTRY_MEAT_HPP
#define EPIWORLD_REGISTRY_MEAT_HPP

inline ParamInfo ParamInfo::probability(
    std::string arg, std::string name, double value, std::string description
)
{
    return {std::move(arg), std::move(name), value, 0.0, 1.0, false,
        std::move(description)};
}

inline ParamInfo ParamInfo::duration(
    std::string arg, std::string name, double value, std::string description
)
{
    return {std::move(arg), std::move(name), value, 0.0,
        std::numeric_limits< double >::infinity(), false,
        std::move(description)};
}

inline ParamInfo ParamInfo::period(
    std::string arg, std::string name, double value, std::string description
)
{
    return {std::move(arg), std::move(name), value, -1.0,
        std::numeric_limits< double >::infinity(), true,
        std::move(description)};
}

inline ParamInfo ParamInfo::count(
    std::string arg, std::string name, double value, std::string description
)
{
    return {std::move(arg), std::move(name), value, 0.0,
        static_cast< double >(std::numeric_limits< epiworld_fast_uint >::max()),
        true, std::move(description)};
}

inline ParamInfo ParamInfo::agents(double value)
{
    ParamInfo p = count("n", "", value, "Number of agents.");
    p.min = 1.0;
    return p;
}

inline ParamInfo ParamInfo::prevalence(double value)
{
    return probability("prevalence", "", value,
        "Proportion of agents infected at the start.");
}

inline ModelArgs::ModelArgs(
    std::map< std::string, double > values_, std::string vname_
) : values(std::move(values_)), vname(std::move(vname_)) {}

inline double ModelArgs::operator()(const std::string & arg) const
{
    auto it = values.find(arg);
    if (it == values.end())
        throw std::logic_error(
            "The model's builder reads \"" + arg +
            "\", which is not one of its arguments."
        );
    return it->second;
}

inline const std::string & ModelArgs::virus_name() const noexcept
{
    return vname;
}

inline ModelRegistry & ModelRegistry::add(ModelInfo info)
{

    if (info.id.empty())
        throw std::invalid_argument("A registered model needs an id.");

    for (const auto & m : models)
        if (m.id == info.id)
            throw std::invalid_argument(
                "The registry already has a model \"" + info.id + "\"."
            );

    if (!info.build)
        throw std::invalid_argument(
            "The model \"" + info.id + "\" has no builder."
        );

    for (size_t i = 0u; i < info.params.size(); ++i)
    {
        const auto & p = info.params[i];

        for (size_t j = 0u; j < i; ++j)
            if (info.params[j].arg == p.arg)
                throw std::invalid_argument(
                    "The model \"" + info.id + "\" lists the argument \"" +
                    p.arg + "\" twice."
                );

        if ((p.value < p.min) || (p.value > p.max))
            throw std::invalid_argument(
                "The default of \"" + p.arg + "\" in the model \"" + info.id +
                "\" is outside its limits."
            );
    }

    models.push_back(std::move(info));
    return *this;

}

inline const std::vector< ModelInfo > &
ModelRegistry::list_models() const noexcept
{
    return models;
}

inline const ModelInfo & ModelRegistry::get_model_info(
    const std::string & id
) const
{

    for (const auto & m : models)
        if (m.id == id)
            return m;

    std::string known;
    for (const auto & m : models)
        known += (known.empty() ? "" : ", ") + m.id;

    throw std::invalid_argument(
        "Unknown model \"" + id + "\". The registry has: " + known + "."
    );

}

inline std::unique_ptr< Model > ModelRegistry::make_model(
    const std::string & id,
    const std::map< std::string, double > & values,
    const std::string & virus_name
) const
{

    const auto & info = get_model_info(id);

    // Every value given must be an argument of the model
    for (const auto & kv : values)
    {
        bool found = false;
        for (const auto & p : info.params)
            if (p.arg == kv.first)
            {
                found = true;
                break;
            }

        if (!found)
        {
            std::string known;
            for (const auto & p : info.params)
                known += (known.empty() ? "" : ", ") + p.arg;

            throw std::invalid_argument(
                "The model \"" + id + "\" has no argument \"" + kv.first +
                "\". Its arguments are: " + known + "."
            );
        }
    }

    // Complete with the defaults, and check
    std::map< std::string, double > complete;
    for (const auto & p : info.params)
    {
        auto it = values.find(p.arg);
        const double v = (it == values.end()) ? p.value : it->second;

        if (!(v >= p.min && v <= p.max))
            throw std::invalid_argument(
                "The argument \"" + p.arg + "\" of the model \"" + id +
                "\" must be in [" + std::to_string(p.min) + ", " +
                std::to_string(p.max) + "]; got " + std::to_string(v) + "."
            );

        if (p.integer && (v != std::floor(v)))
            throw std::invalid_argument(
                "The argument \"" + p.arg + "\" of the model \"" + id +
                "\" must be a whole number; got " + std::to_string(v) + "."
            );

        complete[p.arg] = v;
    }

    if (!virus_name.empty() && info.virus_name.empty())
        throw std::invalid_argument(
            "The model \"" + id + "\" names its own virus; it takes no "
            "virus name."
        );

    return info.build(ModelArgs(
        std::move(complete),
        virus_name.empty() ? info.virus_name : virus_name
    ));

}

#endif
