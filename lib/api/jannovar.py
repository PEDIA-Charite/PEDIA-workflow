'''Client Interface to Jannovar VCF converter server.
'''
import socket
import select
import logging
import io
import typing

from contextlib import contextmanager

from lib.vcf_jannovar import jannovar_vcf_to_table
from lib.singleton import LazyConfigure

LOGGER = logging.getLogger(__name__)


class JannovarClient(LazyConfigure):
    '''Implements a basic stream socket client.
    Built after simple toy implementation here:
    https://docs.python.org/3/howto/sockets.html
    '''

    def __init__(
            self,
    ):
        super().__init__()
        self.url = None
        self.port = None

    def configure(
            self,
            url: str = "localhost",
            port: int = 8888,
    ):
        super().configure()
        self.url = url
        self.port = port

    def create_vcf(
            self,
            variants: [str],
            zygosity: str,
            case_id: str,
    ) -> typing.Union["pandas.DataFrame", str]:
        '''Create pandas dataframe with vcf information.'''
        status, vcf_text = self.process_variants(variants)
        # return error
        if status < 0:
            return vcf_text

        with io.StringIO(vcf_text) as reader:
            vcf_table = jannovar_vcf_to_table(
                reader, case_id, zygosity, variants
            )

        return vcf_table

    def process_variants(self, variants: [str]) -> (int, str):
        '''Submit hgvs vcf file from server.'''
        msg = "\n".join(variants) + "\n"

        msg = "{}\n{}".format(len(msg), msg)

        byte_msg = msg.encode("utf-8")

        with self._connect() as con:
            try:
                self._send(con, byte_msg)
                ready = select.select([con], [], [], 20)
                if ready[0]:
                    status, raw_data = self._recv(con)
                    data = raw_data.decode("utf-8")
            except RuntimeError as error:
                LOGGER.warning(error)
                data = None
        return status, data

    def can_connect(self):
        '''Test whether server can be reached.'''
        try:
            with self._connect():
                pass
        except ConnectionRefusedError:
            return False
        return True

    @contextmanager
    def _connect(self):
        '''Handle connection to jannovar server.'''
        connection = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        connection.connect((self.url, self.port))
        yield connection
        connection.close()

    @staticmethod
    def _send(sock: socket.socket, msg: bytes):
        sent_bytes = 0

        while sent_bytes < len(msg):
            sent = sock.send(msg[sent_bytes:])
            if sent == 0:
                raise RuntimeError("connection broken")
            sent_bytes += sent

    @staticmethod
    def _recv(sock: socket.socket) -> (int, bytes):
        # infer size from first line
        initial = sock.recv(2048)
        msglen, status, chunk = initial.split(b"\n", 2)
        msglen = int(msglen.decode("utf-8"))
        status = int(status.decode("utf-8"))
        bytes_recd = len(chunk)
        chunks = [chunk]

        while bytes_recd < msglen:
            chunk = sock.recv(min(msglen - bytes_recd, 2048))

            if chunk == b"":
                raise RuntimeError("connection broken")

            chunks.append(chunk)
            bytes_recd += len(chunk)

        return status, b"".join(chunks)
