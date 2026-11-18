#ifndef EPIWORLD_VIRUS_HPP
#define EPIWORLD_VIRUS_HPP

class Agent;

class Virus;

class Model;

class Tool;

/**
 * @brief Virus
 * 
 * @details
 * A virus is identified by its id, which the model assigns: `Model::add_virus()`
 * records a new variant (the founder of a lineage), and so does every
 * mutation (a child of the variant it changed). Tools can reduce what a
 * virus does, so its effects are a function of the agent.
 */
class Virus {
    friend class Agent;
    friend class Model;
    friend class DataBase;
    friend class Tool;
    friend class ContactLayer;
    friend struct DiseaseTargets;
private:
    
    Agent * agent = nullptr;


    std::string virus_name = "unknown virus";
    int date = -99;
    int id   = -99;    
    int lineage_id = -99; ///< Id of the founding virus (kept across mutations).
    int disease = -1;     ///< The model's disease it belongs to (see `Model::add_virus()`)
    uint64_t disease_bit = 0u; ///< Its bit in a target mask (none: not in a model)
    MutFun          mutation                 = nullptr;

    /// Parameters this variant reads in place of others (see `use_param()`)
    std::vector< std::pair< ParamRef, ParamRef > > param_redirects;

    /// Its natural history (shared by its copies; copied on the first edit)
    std::shared_ptr< NaturalHistory > nh_def = nullptr;
    NaturalHistory & nh_edit();
        
public:

    #ifdef EPI_DEBUG_VIRUS
    static std::atomic<int> counter_construct;          // Default and parameterized constructors
    static std::atomic<int> counter_copy_construct;     // Copy constructor
    static std::atomic<int> counter_move_construct;     // Move constructor
    static std::atomic<int> counter_copy_assign;        // Copy assignment
    static std::atomic<int> counter_move_assign;        // Move assignment
    static std::atomic<int> counter_destruct;           // Destructor
    #endif

    Virus();

    Virus(std::string name = "unknown virus");
    
    
    #ifdef EPI_DEBUG_VIRUS
    
    // Copy and move operations for debugging
    Virus(const Virus& other);                    // Copy constructor
    Virus(Virus&& other) noexcept;                // Move constructor
    Virus& operator=(const Virus& other);   // Copy assignment
    Virus& operator=(Virus&& other) noexcept; // Move assignment

    ~Virus();
    #endif
    
    /**
     * @brief Mutation
     * @details The function returns `true` when it changed the virus (it may
     * change its name or, with `use_param()`, the parameters it reads). The
     * model then records the virus as a new variant: a new id, whose parent is
     * the id it had, in the same lineage. Each mutation is a new variant, even
     * if two of them made the same change.
     */
    ///@{
    void mutate(Model * model);
    void set_mutation(MutFun fun);
    ///@}
    
    
    Agent * get_agent();
    void set_agent(Agent * p);
    
    void set_date(int d);
    int get_date() const;

    void set_id(int idx);
    int get_id() const;

    /**
     * @brief Id of the virus lineage
     * 
     * The id assigned to the founding virus when it was added to the model
     * (see `Model::add_virus()`). Mutations change the virus id but keep
     * the lineage id, so all variants of a virus share it.
     * 
     * @return int The lineage id (-99 if the virus has not been recorded).
     */
    int get_lineage_id() const;

    void set_name(std::string name);
    std::string get_name() const;

    bool operator==(const Virus & other) const;
    bool operator!=(const Virus & other) const {return !operator==(other);};

    void print() const;


    /**
     * @name Natural history
     * @brief The disease's stages and the transitions between them, declared
     * on the virus (see `NaturalHistory` for each method).
     * @details A natural history is a disease. A virus built from scratch
     * (with `add_stage()` and the others) starts a new disease; a copy of it
     * (a strain, see `strain()`, or a variant from a mutation) belongs to the
     * same disease: it keeps the stages and transitions, and can read other
     * parameters with `use_param()` (see `Model::add_virus()`).
     */
    ///@{
    int add_stage(std::string label, unsigned int flags = 0u);
    void set_contact(int s, epiworld_double weight);
    TransitionRef infection(Nodes from, int to, Value per_contact);
    TransitionRef rate(Nodes from, int to, Value daily_probability);
    TransitionRef after(int from, int to, Duration d);
    TransitionRef external(Nodes from, int to, std::string label);

    /// The natural history (`nullptr` if the virus has none).
    const NaturalHistory * get_natural_history() const { return nh_def.get(); }

    /**
     * @brief A strain of this virus: a copy named `name`, of the same disease.
     * @details It shares the natural history, so `Model::add_virus()` adds it
     * to the disease of this virus, and can read other parameters
     * (`use_param()`). It is not in a model until added.
     */
    Virus strain(std::string name) const;

    /// The index of its disease in its model (`Model::disease_id()`; -1 if it
    /// is not in a model).
    int get_disease() const { return disease; }
    ///@}

    /**
     * @brief This variant reads the model parameter `replacement` wherever a
     * transition of the natural history reads `param`.
     * @details Every parameter lives in the model; a variant changes a value
     * by pointing to another parameter (e.g., `"Recovery rate (Delta)"`).
     */
    void use_param(std::string param, std::string replacement);

    /// The parameter this variant reads in place of `id`.
    ParamId param_redirect(ParamId id, const Model & model) const;

    /// Whether it reads any parameter in place of another (`use_param()`).
    bool has_param_redirects() const { return !param_redirects.empty(); }

    /// Throws if a parameter it reads (see `use_param()`) is not in `model`.
    void check_params(const Model & model) const;

    virtual ~Virus() = default;
    virtual std::unique_ptr<Virus> clone_ptr() const;

};

#endif