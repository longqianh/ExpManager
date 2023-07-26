function tcp_client=create_tcp_client(address,port,callback_fun)
    % create server in C++ first
    read_count=1;
    tcp_client=tcpclient(address,port);
    configureCallback(tcp_client,"byte",read_count,callback_fun);
end