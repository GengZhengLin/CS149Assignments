// ChatServer

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChatServer {
    private static final Charset utf8 = Charset.forName("UTF-8");

    private static final String OK = "200 OK";
    private static final String NOT_FOUND = "404 NOT FOUND";
    private static final String HTML = "text/html";
    private static final String TEXT = "text/plain";
    private static final String DEFAULT_ROOM = "DEFAULT";

    private static final Pattern PAGE_REQUEST
        = Pattern.compile("GET /([\\p{Alnum}]*/?) HTTP.*");
    private static final Pattern PULL_REQUEST
        = Pattern.compile("POST /([\\p{Alnum}]*)/?pull\\?last=([0-9]+) HTTP.*");
    private static final Pattern PUSH_REQUEST
        = Pattern.compile("POST /([\\p{Alnum}]*)/?push\\?msg=([^ ]*) HTTP.*");

    private static final String CHAT_HTML;
    static {
        try {
            CHAT_HTML = getFileAsString("../index.html");
        } catch (final IOException xx) {
            throw new Error("unable to start server", xx);
        }
    }

    private final int port;
    private final Map<String,ChatState> stateByName
        = new HashMap<String,ChatState>();
    private ChatState chatStateAll;
    
    class ServerThread extends Thread {
    	Socket connection;
    	public void setConnection(final Socket connection){
    		this.connection = connection;
    	}
    	
    	public void run(){
    		try{
    			ChatServer.this.handle(connection);
    		}
    		catch(Exception e){
    			System.err.println("Caught IOException: " + e.getMessage());
    			System.exit(-1);
    		} 
    		finally{
    			synchronized(ChatServer.this){
    				ChatServer.this.notify();
    			}
    		}
    	}
    }
    
    private ServerThread[] serverThreads = new ServerThread[8];

    /**
     * Constructs a new {@link ChatServer} that will service requests
     * on the specified <code>port</code>. <code>state</code> will be
     * used to hold the current state of the chat.
     */
    public ChatServer(final int port) throws IOException {
        this.port = port;
        chatStateAll = new ChatState("all");
        stateByName.put("all", chatStateAll);
    }

    public void runForever() throws IOException {
        @SuppressWarnings("resource")
		final ServerSocket server = new ServerSocket(port);
        while (true) {
        	final Socket connection = server.accept();
        	boolean getJobWorking = false;
        	while (!getJobWorking){
                for(int i = 0; i < serverThreads.length; i++){
                	if (serverThreads[i] == null || !serverThreads[i].isAlive()){
                		serverThreads[i] = new ServerThread();
                		serverThreads[i].setConnection(connection);
                		serverThreads[i].start();
                		getJobWorking = true;
                		break;
                	}
                }
                if (!getJobWorking){
                	try {
                		synchronized(this){
                			wait();
                		}
                	}
                	catch(InterruptedException e) {
                		e.printStackTrace();
                	}
                }
        	}
        }
    }
    
    private static String replaceEmptyWithDefaultRoom(final String room) {
    	if (room.isEmpty()) {
    		return DEFAULT_ROOM;
    	}
    	return room;
    }

    private void handle(final Socket connection) throws IOException, InterruptedException {
        try {
            final BufferedReader xi
                = new BufferedReader(new InputStreamReader(connection.getInputStream()));
            final OutputStream xo = new BufferedOutputStream(connection.getOutputStream());

            final String request = xi.readLine();
            System.out.println(Thread.currentThread() + ": " + request);

            Matcher m;
            if (request == null) {
                sendResponse(xo, NOT_FOUND, TEXT, "Empty request.");
            } else if (PAGE_REQUEST.matcher(request).matches()) {
                sendResponse(xo, OK, HTML, CHAT_HTML);
            } else if ((m = PULL_REQUEST.matcher(request)).matches()) {
                String room = replaceEmptyWithDefaultRoom(m.group(1));
                final long last = Long.valueOf(m.group(2));
                sendResponse(xo, OK, TEXT, getState(room).recentMessages(last));
            } else if ((m = PUSH_REQUEST.matcher(request)).matches()) {
                String room = replaceEmptyWithDefaultRoom(m.group(1));
                final String msg = m.group(2);
                if (room.equals("all")){
                	Collection<ChatState> rooms = null;
                	synchronized(stateByName){
                		rooms = stateByName.values();               	
	                	for (ChatState value : rooms){
	                		synchronized(value){
	                			value.addMessage(msg);
	                		}
	                	}
                	}
                }
                else {
                	ChatState chatRoom = getState(room);
                	ChatState chatAll = getState("all");
                	Thread currentThread = Thread.currentThread();
            		synchronized(chatRoom){
                    	chatRoom.addMessage(msg);
                    	chatRoom.addingMsgThreads.add(currentThread);
            		}
            		Thread allRoomThread = new Thread(new Runnable(){
            			public void run(){
        					try {
                    			synchronized (chatRoom){
                    				if (chatRoom.addingMsgThreads.peek() != currentThread){
    									chatRoom.wait();
                    				}
                        			synchronized (chatAll) {              				
                        				chatAll.addMessage(msg);
            						}
                    				chatRoom.addingMsgThreads.remove();
                    				chatRoom.notifyAll();
                    			}
							} catch (Exception e) {
								e.printStackTrace();
							}
            			}
            		});
            		allRoomThread.start();
            		sendResponse(xo, OK, TEXT, "ack");
                }                
            } else {
                sendResponse(xo, NOT_FOUND, TEXT, "Malformed request.");
            }
        } finally {
        	connection.close();
        }
    }
    
    /**
     * Writes a minimal but valid HTTP response to
     * <code>output</code>.
     */
    private static void sendResponse(final OutputStream xo,
                                     final String status,
                                     final String contentType,
                                     final String content) throws IOException {
        final byte[] data = content.getBytes(utf8);
        final String headers =
            "HTTP/1.0 " + status + "\r\n" +
            "Content-Type: " + contentType + "; charset=utf-8\r\n" +
            "Content-Length: " + data.length + "\r\n\r\n";

        xo.write(headers.getBytes(utf8));
        xo.write(data);
        xo.flush();

        System.out.println(Thread.currentThread() + ": replied with " + data.length + " bytes");
    }

    private ChatState getState(final String room) {
    	synchronized(stateByName){
            ChatState state = stateByName.get(room);
            if (state == null) {
                state = new ChatState(room);
                stateByName.put(room, state);
            }
            return state;
    	}       
    }

    /**
     * Reads the resource with the specified path as a string, and
     * then returns the string.
     */
    private static String getFileAsString(final String path)
        throws IOException {
    	byte[] fileBytes = Files.readAllBytes(Paths.get(path));
    	return new String(fileBytes, utf8);
    }

    /**
     * Runs a chat server, with a default port of 8080.
     */
    public static void main(final String[] args) throws IOException {
        final int port = args.length == 0 ? 8080 : Integer.parseInt(args[0]);
        new ChatServer(port).runForever();
    }
}
