#include "colvarscript.h"

#include "colvarmodule.h"
#include "colvarproxy.h"

#include <cstring>
#include <sstream>

namespace {

template <typename T>
std::string value_to_string(const T& value)
{
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

template <typename T>
std::string vector_to_string(const std::vector<T>& values)
{
    std::ostringstream stream;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            stream << ' ';
        }
        stream << values[i];
    }
    return stream.str();
}

} // namespace

colvarscript::colvarscript(colvarproxy* proxy, colvarmodule* module)
    : proxy_(proxy),
      cvmodule(module),
      cmd_names(nullptr)
{
}

colvarscript::~colvarscript() = default;

int colvarscript::run(int, unsigned char* const*)
{
    return unsupported_op();
}

int colvarscript::set_result_str(std::string const& value)
{
    str_result_ = value;
    return COLVARSCRIPT_OK;
}

int colvarscript::clear_str_result()
{
    str_result_.clear();
    return COLVARSCRIPT_OK;
}

void colvarscript::add_error_msg(std::string const& message)
{
    proxy_->add_error_msg(message);
}

std::string colvarscript::get_cmd_prefix(Object_type)
{
    return "cv";
}

char const* colvarscript::get_command_help(char const*)
{
    return "Scripting is not supported by the OpenMM Colvars interface.";
}

char const* colvarscript::get_command_rethelp(char const*)
{
    return "";
}

char const* colvarscript::get_command_arghelp(char const*, int)
{
    return "";
}

int colvarscript::get_command_n_args_min(char const*)
{
    return 0;
}

int colvarscript::get_command_n_args_max(char const*)
{
    return 0;
}

char const* colvarscript::get_command_full_help(char const*)
{
    return "Scripting is not supported by the OpenMM Colvars interface.";
}

std::string colvarscript::get_cmdline_help_summary(Object_type)
{
    return "Scripting is not supported by the OpenMM Colvars interface.";
}

std::string colvarscript::get_command_cmdline_syntax(Object_type, command)
{
    return "";
}

std::string colvarscript::get_command_cmdline_help(Object_type, std::string const&)
{
    return "Scripting is not supported by the OpenMM Colvars interface.";
}

int colvarscript::unsupported_op()
{
    add_error_msg("Error: scripting is not supported by the OpenMM Colvars interface.\n");
    return COLVARSCRIPT_ERROR;
}

char* colvarscript::obj_to_str(unsigned char* obj)
{
    return reinterpret_cast<char*>(obj);
}

std::vector<std::string> colvarscript::obj_to_str_vector(unsigned char*)
{
    return {};
}

int colvarscript::set_result_int(int const& value, unsigned char* obj)
{
    return set_result_text_from_str(value_to_string(value), obj);
}

int colvarscript::set_result_int_vec(std::vector<int> const& value, unsigned char* obj)
{
    return set_result_text_from_str(vector_to_string(value), obj);
}

int colvarscript::set_result_long_int(long int const& value, unsigned char* obj)
{
    return set_result_text_from_str(value_to_string(value), obj);
}

int colvarscript::set_result_long_int_vec(std::vector<long int> const& value, unsigned char* obj)
{
    return set_result_text_from_str(vector_to_string(value), obj);
}

int colvarscript::set_result_real(cvm::real const& value, unsigned char* obj)
{
    return set_result_text_from_str(cvm::to_str(value), obj);
}

int colvarscript::set_result_rvector(cvm::rvector const& value, unsigned char* obj)
{
    return set_result_text_from_str(cvm::to_str(value), obj);
}

int colvarscript::set_result_colvarvalue(colvarvalue const& value, unsigned char* obj)
{
    return set_result_text_from_str(value.to_simple_string(), obj);
}

int colvarscript::set_result_colvarvalue_vec(std::vector<colvarvalue> const& value,
                                             unsigned char* obj)
{
    std::ostringstream stream;
    for (size_t i = 0; i < value.size(); ++i) {
        if (i > 0) {
            stream << ' ';
        }
        stream << value[i].to_simple_string();
    }
    return set_result_text_from_str(stream.str(), obj);
}

int colvarscript::init_commands()
{
    return COLVARSCRIPT_OK;
}

int colvarscript::init_command(colvarscript::command const&,
                               char const*,
                               char const*,
                               int,
                               int,
                               char const*,
                               int (*)(void*, int, unsigned char* const*))
{
    return COLVARSCRIPT_OK;
}

int colvarscript::proc_features(colvardeps*, int, unsigned char* const*)
{
    return unsupported_op();
}

int colvarscript::set_result_text_from_str(std::string const& value, unsigned char* obj)
{
    if (obj != nullptr) {
        std::strcpy(reinterpret_cast<char*>(obj), value.c_str());
        return COLVARSCRIPT_OK;
    }
    str_result_ = value;
    return COLVARSCRIPT_OK;
}

extern "C" int run_colvarscript_command(colvarscript* script,
                                        int objc,
                                        unsigned char* const objv[])
{
    if (script == nullptr) {
        return COLVARSCRIPT_ERROR;
    }
    return script->run(objc, objv);
}

extern "C" const char* get_colvarscript_result()
{
    return "";
}
