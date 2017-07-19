#ifndef SETTINGS_H
#define SETTINGS_H

struct Settings {
    unsigned redundancy_check_skip_recursions{};
    bool allow_weight_one{};
    bool allow_weight_one_half{};
    bool allow_weights_sum_one{};
    bool print_last_recursion_statistics{};
    bool print_candidates{};
    bool debug_ignore_symmetries{};
    bool debug_disable_lex_order{};
    bool ip_check{};
    bool count_weight_systems{};
};

extern Settings g_settings;

#endif
