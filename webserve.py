import time
import SimpleHTTPServer
import SocketServer
HOST_NAME = "localhost"
PORT = 8080

Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

httpd = SocketServer.TCPServer(("", PORT), Handler)

try:
	print "serving at port ", PORT
	httpd.serve_forever()
except KeyboardInterrupt:
	pass
httpd.server_close()
print time.asctime(), "Server Stops - %s:%s" % (HOST_NAME, PORT)
