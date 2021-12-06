function output = encode_binary(binary_array)
output = num2str(binary_array);
output(isspace(output)) = '';
end